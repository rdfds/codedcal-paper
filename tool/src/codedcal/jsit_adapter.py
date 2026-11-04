"""Disclosed Python adapter of JSIT's published MATLAB SGL/FISTA kernel.

Upstream 3aa4e398760f4b576cd592eea46a37bdeef23073 is MIT licensed.
This offers an exact published-kernel replay and a separately labelled
standard-FISTA correction. It is not the full native pipeline: image
preparation, patching and filtering are explicit.
"""
import numpy as np
from scipy import sparse
from scipy.ndimage import label, center_of_mass


def psf_matrix(side,scale=3,sigma=1.25):
    """Exact separable form of getPsfMat2's finite normalized fspecial kernel.

    MATLAB's even Gaussian kernel has half-pixel centers; detected pixels
    sum scale-by-scale high-resolution cells with zero outside the patch.
    """
    if side<1 or int(side)!=side or int(scale)!=scale or scale<1 or sigma<=0:
        raise ValueError('positive integer dimensions and positive sigma required')
    n=int(side*scale);axis=np.arange(2*n)-n+.5
    kernel=np.exp(-axis**2/(2*sigma*sigma));kernel/=kernel.sum()
    # Rows are low-resolution detectors; columns are high-resolution sources.
    fine=np.arange(n)[:,None]-np.arange(n)[None,:]+n
    one=kernel[fine].reshape(side,scale,n).sum(axis=1)
    one[np.abs(one)<1e-14]=0
    return sparse.kron(sparse.csr_matrix(one),sparse.csr_matrix(one),format='csr')


def prox_sgl(x,threshold,alpha=.5):
    if threshold<0 or not 0<=alpha<=1:raise ValueError('invalid proximal parameters')
    soft=np.sign(x)*np.maximum(np.abs(x)-alpha*threshold,0.)
    norms=np.linalg.norm(soft,axis=1)
    shrink=np.maximum(1-(1-alpha)*threshold/np.maximum(norms,np.finfo(float).tiny),0)
    return soft*shrink[:,None]


def fista(A,C,Y,penalty,iterations=10,alpha=.5,eK=None,eM=None,variant="standard"):
    """SGL model with standard FISTA, or exact published kernel replay.

    The published kernel fails a zero-penalty identity-operator oracle because
    its step starts at X with a gradient at z and uses t/t_previous momentum.
    `published` preserves that behavior for native equivalence only.
    `standard` explicitly corrects both formulas; it is an adapter, not a
    claim to run the unmodified full native JSIT production pipeline.
    """
    if variant not in ["standard","published"]:raise ValueError("unknown FISTA variant")
    C=np.asarray(C,float);Y=np.asarray(Y,float)
    if C.ndim!=2 or Y.shape!=(A.shape[0],C.shape[1]) or not np.isfinite(Y).all():
        raise ValueError('finite matching measurements required')
    if penalty<0 or iterations<1:raise ValueError('invalid optimization parameters')
    if eK is None:eK=float(np.linalg.eigvalsh(C.T@C)[-1])
    if eM is None:eM=float(np.linalg.eigvalsh((A@A.T).toarray() if sparse.issparse(A) else A@A.T)[-1])
    L=eK*eM
    if not np.isfinite(L) or L<=0:raise ValueError('positive finite Lipschitz constant required')
    X=np.zeros((A.shape[1],C.shape[0]));z=X.copy();t=1.;losses=[]
    for _ in range(iterations):
        residual=A@(z@C)-Y
        gradient=A.T@residual@C.T
        previous=X;X=prox_sgl((X if variant=="published" else z)-gradient/L,penalty/L,alpha)
        old_t=t;t=.5*(1+np.sqrt(1+4*t*t))
        momentum=t/old_t if variant=="published" else (old_t-1)/t
        z=X+momentum*(X-previous)
        error=A@(X@C)-Y
        losses.append(float(.5*np.sum(error*error)+penalty*(alpha*np.abs(X).sum()+(1-alpha)*np.linalg.norm(X,axis=1).sum())))
    if not np.isfinite(X).all():raise ValueError('nonfinite JSIT adapter solution')
    return X,dict(iterations=iterations,penalty=penalty,alpha=alpha,objective=losses,
        variant=variant,extrapolation="published_t_over_previous_t" if variant=="published" else "standard_previous_t_minus_one_over_t",signed_coefficients=True)


def component_calls(values,threshold=0.,min_area=2,scale=1,border=8):
    """One maximal-absolute word per location, then 4-connected components.

    JSIT hard threshold is signed (negative values below threshold are removed).
    Thus nonnegative thresholds exclude negative coefficients before winning.
    Winner restriction loses exactly co-located different genes; report this
    limitation rather than claiming an unrestricted source demixer.
    """
    values=np.asarray(values,float)
    if values.ndim!=3 or not np.isfinite(values).all() or threshold<0 or min_area<1:
        raise ValueError('finite density, nonnegative threshold and positive area required')
    v=np.where(values>=threshold,values,0)
    winners=np.argmax(np.abs(v),axis=-1);amplitude=np.take_along_axis(v,winners[...,None],axis=-1)[...,0]
    xy=[];words=[];areas=[];strength=[]
    h,w=winners.shape
    for k in np.unique(winners[amplitude>0]):
        components,n=label((winners==k)&(amplitude>0))
        size=np.bincount(components.ravel())
        total=np.bincount(components.ravel(),weights=amplitude.ravel())
        ids=np.flatnonzero(size>=min_area);ids=ids[ids>0]
        for i,(y,x) in zip(ids,center_of_mass(np.ones_like(amplitude),components,ids)):
            # Native centroid uses 1-based coordinates, divides by scale.
            # Convert that convention back to 0-based camera coordinates.
            x=(x+1)/scale-1;y=(y+1)/scale-1
            if border<=x<w/scale-border and border<=y<h/scale-border:
                xy.append([x,y]);words.append(k);areas.append(int(size[i]));strength.append(float(total[i]/size[i]))
    return np.asarray(xy).reshape(-1,2),np.asarray(words,int),dict(areas=areas,mean_coefficient=strength)


def decode_density(stack,words,penalty=.5,sigma=3.,scale=3,patch_side=16,iterations=10,variant="standard"):
    """Independent patches with an explicit resolution/PSF adapter."""
    stack=np.asarray(stack,float);r,h,w=stack.shape
    if h%patch_side or w%patch_side:raise ValueError('patches must tile exactly')
    A=psf_matrix(patch_side,scale,sigma);C=np.asarray(words,float)
    eK=float(np.linalg.eigvalsh(C.T@C)[-1]);eM=float(np.linalg.eigvalsh((A@A.T).toarray())[-1])
    result=np.zeros((h*scale,w*scale,len(words)),np.float32);diags=[]
    for y in range(0,h,patch_side):
        for x in range(0,w,patch_side):
            Y=np.moveaxis(stack[:,y:y+patch_side,x:x+patch_side],0,-1).reshape(-1,r,order='F')
            X,diag=fista(A,C,Y,penalty,iterations,eK=eK,eM=eM,variant=variant)
            result[y*scale:(y+patch_side)*scale,x*scale:(x+patch_side)*scale]=X.reshape(patch_side*scale,patch_side*scale,len(words),order='F')
            diags.append(diag)
    return result,dict(patches=len(diags),first_patch=diags[0],last_patch=diags[-1],scale=scale,
        high_resolution_sigma=sigma,patch_side=patch_side,upstream_commit='3aa4e398760f4b576cd592eea46a37bdeef23073',
        native_full_workflow=False,preprocessing_not_part_of_this_function=True)
