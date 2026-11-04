"""Nonnegative sparse decoding with a preprocessing-matched forward operator.

The Gaussian-difference preparation is linear and must also act on the PSF.
A pure Gaussian model of these signed images instead forces negative halos
into latent source coefficients. This is a disclosed model adapter, not native
JSIT. Pixel/code weights must be trained outside the scored regions.
"""
import numpy as np
from scipy.ndimage import gaussian_filter
from .jsit_adapter import psf_matrix


def effective_psf_matrix(side,scale=3,raw_sigma=3.,foreground_sigma=.5,background_sigma=3.):
    base=psf_matrix(side,scale,raw_sigma).toarray()
    images=base.reshape(side,side,base.shape[1],order='F')
    prepared=gaussian_filter(images,(foreground_sigma,foreground_sigma,0))-gaussian_filter(images,(background_sigma,background_sigma,0))
    return prepared.reshape(base.shape,order='F')


def prox_nonnegative_sgl(x,threshold,alpha=.5):
    if threshold<0 or not 0<=alpha<=1:raise ValueError('invalid penalty')
    soft=np.maximum(x-alpha*threshold,0)
    norm=np.linalg.norm(soft,axis=1)
    return soft*np.maximum(1-(1-alpha)*threshold/np.maximum(norm,np.finfo(float).tiny),0)[:,None]


def solve(A,C,Y,penalty=.1,iterations=150,alpha=.5):
    A=np.asarray(A,float);C=np.asarray(C,float);Y=np.asarray(Y,float)
    if Y.shape!=(A.shape[0],C.shape[1]) or not all(np.isfinite(x).all() for x in [A,C,Y]) or np.any(C<0):raise ValueError('invalid finite model dimensions')
    L=float(np.linalg.eigvalsh(A@A.T)[-1]*np.linalg.eigvalsh(C.T@C)[-1])
    if L<=0 or penalty<0 or iterations<1:raise ValueError('invalid optimization state')
    x=np.zeros((A.shape[1],len(C)));z=x.copy();t=1.;objectives=[]
    for _ in range(iterations):
        gradient=A.T@(A@(z@C)-Y)@C.T
        previous=x;x=prox_nonnegative_sgl(z-gradient/L,penalty/L,alpha)
        old=t;t=(1+np.sqrt(1+4*t*t))/2;z=x+(old-1)/t*(x-previous)
        residual=A@(x@C)-Y
        objectives.append(float(.5*np.sum(residual*residual)+penalty*(alpha*x.sum()+(1-alpha)*np.linalg.norm(x,axis=1).sum())))
    gradient=A.T@(A@(x@C)-Y)@C.T
    fixed=prox_nonnegative_sgl(x-gradient/L,penalty/L,alpha)
    return x,dict(iterations=iterations,penalty=penalty,objective=objectives,proximal_fixed_point_relative=float(np.linalg.norm(x-fixed)/max(np.linalg.norm(x),1e-12)),
       nonnegative_sources=True,forward_operator='Gaussian emitter PSF followed by same Gaussian-difference preparation',lipschitz=L)
