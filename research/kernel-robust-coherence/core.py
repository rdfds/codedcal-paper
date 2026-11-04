"""Research-only matched-four-count coherence. No physical validity certificate."""
import itertools,math
import numpy as np

COEFF={1:1.,2:-1/3,3:1/3,4:-1.}
NORM={2:1/3,3:1/6,4:5/24}
LABELS=np.array(list(itertools.product(range(4),repeat=4)))
KERNEL=np.array([COEFF[len(set(c))] for c in LABELS]).reshape((4,)*4)


def validate(z,n):
    z=np.asarray(z,complex);n=np.asarray(n)
    if z.ndim<1 or z.shape!=n.shape or z.shape[-1]!=4 or np.iscomplexobj(n):raise ValueError('matching real-count/complex-phase arrays ending in four panels required')
    n=n.astype(float)
    if not np.isfinite(z).all() or not np.isfinite(n).all() or np.any(n<0):raise ValueError('finite nonnegative counts required')
    if np.any(abs(z)>n+1e-10*np.maximum(1,n)):raise ValueError('phase sum exceeds count capacity')
    return z,n


def validate_codes(codes,sources):
    codes=np.asarray(codes)
    if codes.shape!=(sources,) or not np.issubdtype(codes.dtype,np.number) or np.iscomplexobj(codes):raise ValueError('one finite integer code per source required')
    if not np.isfinite(codes).all() or np.any((codes<0)|(codes>3)) or np.any(codes!=np.floor(codes)):raise ValueError('codes must be integers in [0,3]')
    return codes.astype(int)


def collapsed(z,n):
    """Components if the provided phase/count sums have a single source/code."""
    z,n=validate(z,n)
    c=(z[...,0]*n[...,1]+n[...,0]*z[...,1])*(z[...,2].conj()*n[...,3]+n[...,2]*z[...,3].conj())/4
    ra=(z[...,0]*z[...,1].conj()).real*n[...,2]*n[...,3]
    rb=n[...,0]*n[...,1]*(z[...,2]*z[...,3].conj()).real
    return np.stack([c.real,c.imag,ra,rb],axis=-1)


def h4(a,b,c,d):
    return np.einsum('ijkl,...i,...j,...k,...l->...',KERNEL,a,b,c,d,optimize=True)


def coded_components(z,n):
    """Input (...,four codes,four panels); output (...,ReC,ImC,RA,RB)."""
    z,n=validate(z,n)
    if z.shape[-2:]!=(4,4):raise ValueError('four codes and four panels required')
    a,b,c,d=np.moveaxis(z,-1,0);na,nb,nc,nd=np.moveaxis(n,-1,0)
    cr=(h4(a,nb,c.conj(),nd)+h4(a,nb,nc,d.conj())+h4(na,b,c.conj(),nd)+h4(na,b,nc,d.conj()))/4
    return np.stack([cr.real,cr.imag,h4(a,b.conj(),nc,nd).real,h4(na,nb,c,d.conj()).real],axis=-1)


def aggregate(z,n,codes):
    z,n=validate(z,n)
    if z.ndim!=2:raise ValueError('sources by four panels required')
    codes=validate_codes(codes,len(z))
    return tuple(np.stack([v[codes==k].sum(0) for k in range(4)]) for v in (z,n))


def coefficients(z,n,include_count=False):
    """Exact subset coefficients; source-resolved research design input only."""
    z,n=validate(z,n)
    if z.ndim!=2 or len(z)>8:raise ValueError('at most eight resolved sources per research group')
    dim=5 if include_count else 4;coeff={}
    def value(ix):
        zz=z[list(ix)].sum(0);nn=n[list(ix)].sum(0);out=collapsed(zz,nn)
        return np.r_[out,nn.prod()] if include_count else out
    mean=np.zeros(dim)
    for size in range(1,min(4,len(z))+1):
        for subset in itertools.combinations(range(len(z)),size):
            v=value(subset)
            for k in range(1,size):
                for inner in itertools.combinations(subset,k):v-=coeff[inner]
            coeff[subset]=v
            if size==1:mean+=v
    return mean,coeff


def covariance(z,n,include_count=False):
    mean,coeff=coefficients(z,n,include_count);cov=np.zeros((len(mean),len(mean)))
    energies={k:np.zeros_like(cov) for k in NORM}
    for s,b in coeff.items():
        if len(s)>1:energies[len(s)]+=np.outer(b,b)*NORM[len(s)]
    return dict(mean=mean,covariance=sum(energies.values(),np.zeros_like(cov)),energies=energies,coefficients=coeff)


def realize(z,n,codes,include_count=False):
    mean,coeff=coefficients(z,n,include_count);out=mean.copy();codes=validate_codes(codes,len(z))
    for s,b in coeff.items():
        if len(s)>1:out+=b*COEFF[len(set(codes[list(s)]))]
    return out


def collision_aware(z,n,codes):
    """Requires design-known source codes/flag; not inferred from detections."""
    z,n=validate(z,n)
    if z.ndim!=2:raise ValueError('sources by four panels required')
    m=len(z);codes=validate_codes(codes,m)
    if m<1 or m>4:raise ValueError('known occupancy in [1,4] required')
    p=math.prod((4-j)/4 for j in range(m))
    flag=len(set(codes))==m
    return collapsed(z,n).sum(0)/p if flag else np.zeros(4)


def population(n,psi_a,psi_b,m):
    n=np.asarray(n,float);pa=np.asarray(psi_a,complex);pb=np.asarray(psi_b,complex)
    if n.ndim!=2 or n.shape[1]!=4 or pa.shape!=pb.shape or pa.shape!=(len(n),):raise ValueError('incompatible source shapes')
    if not np.isfinite(n).all() or not np.isfinite(pa).all() or not np.isfinite(pb).all() or not np.isfinite(m) or np.any(n<0) or np.any(abs(pa)>1+1e-12) or np.any(abs(pb)>1+1e-12) or abs(m)>1+1e-12:raise ValueError('invalid population inputs')
    w=n.prod(1);c=m*np.sum(w*pa*pb.conj());ra=np.sum(w*abs(pa)**2);rb=np.sum(w*abs(pb)**2)
    return np.array([c.real,c.imag,ra,rb])


def coherence(component):
    x=np.asarray(component,float)
    if x.shape!=(4,) or not np.isfinite(x).all():raise ValueError('four finite real components required')
    re,im,ra,rb=x
    return float(np.hypot(re,im)/np.sqrt(ra*rb)) if ra>0 and rb>0 else None


def witness_vector(g,theta=0.,scale=1.):
    if not 0<=g<=1 or not np.isfinite([g,theta,scale]).all() or scale<=0:raise ValueError('invalid witness')
    return np.array([np.cos(theta),np.sin(theta),-g*scale/2,-g/(2*scale)])


def bound_capacity(n):
    n=np.asarray(n,float)
    if n.ndim!=3 or n.shape[-1]!=4 or np.any(n<0) or not np.isfinite(n).all():raise ValueError('groups,sources,panels counts required')
    return n.sum(1).prod(1)
