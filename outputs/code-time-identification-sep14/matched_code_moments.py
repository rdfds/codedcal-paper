"""Physical four-code diagonal sketches, with no source-assignment inputs.

Input shape (..., 4, 6): per-code [nA,nB,sumXA,sumXB,sumXA2,sumXB2].
This version is one spatial coordinate; add numerators for multiple coordinates.
Raw numerator/denominator are unbiased over independent uniform source codes.
Their ratio is NOT finite-sample unbiased. Do not clip negative estimates.
"""
import itertools
import numpy as np

H = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], float)
PAIRS = tuple(itertools.combinations(range(4),2))
PERMS = tuple(itertools.permutations(range(3)))

def diagonal_moments(coded):
    z = np.einsum('...cf,cb->...fb', coded, H)
    def t2(a,b):
        return np.mean(z[...,a,:]*z[...,b,:],axis=-1)
    def t3(a,b,c):
        v=[z[...,a,:],z[...,b,:],z[...,c,:]]
        return sum(v[0][...,p[0]]*v[1][...,p[1]]*v[2][...,p[2]] for p in PERMS)/6
    def t4(a,b,c,d):
        v=[z[...,a,:],z[...,b,:],z[...,c,:],z[...,d,:]]
        ans=0.
        for x,y in ((0,1),(0,2),(1,2)):
            term=np.prod(np.stack([u[...,x] for u in v]),axis=0)+np.prod(np.stack([u[...,y] for u in v]),axis=0)
            for p in PAIRS:
                term=term-np.prod(np.stack([u[...,x if j in p else y] for j,u in enumerate(v)]),axis=0)
            ans=ans-term/12
        return ans
    return t2,t3,t4

def estimate(coded):
    t2,t3,t4=diagonal_moments(coded)
    a,b,x,y,q,r=range(6)
    num=(t4(b,b,x,x)+t4(a,a,y,y)-2*t4(a,b,x,y)
         -t3(b,x,x)-t3(a,y,y)+2*t3(a,x,y)+2*t3(b,x,y)
         -t3(b,b,q)-t3(a,a,r)+t2(b,q)+t2(a,r)-2*t2(x,y))
    den=t4(a,a,b,b)-t3(a,a,b)-t3(a,b,b)+t2(a,b)
    cross=t2(b,q)+t2(a,r)-2*t2(x,y)
    within_a=t2(a,q)-t2(x,x)
    within_b=t2(b,r)-t2(y,y)
    na,nb=coded[...,0].sum(-1),coded[...,1].sum(-1)
    return dict(numerator=num,denominator=den,cross=cross,cross_den=t2(a,b),
                within_a=within_a,within_b=within_b,
                within_a_den=t2(a,a)-na,within_b_den=t2(b,b)-nb)

def source_polynomial(raw):
    """Audit only: source-level version, shape (..., sources, 6)."""
    a,b,x,y,q,r=np.moveaxis(raw,-1,0)
    cross=b*q+a*r-2*x*y
    wa=a*q-x*x;wb=b*r-y*y
    w=a*(a-1)*b*(b-1)
    d=(a-1)*(b-1)*cross-b*(b-1)*wa-a*(a-1)*wb
    return dict(numerator=d.sum(-1),denominator=w.sum(-1),
                cross=cross.sum(-1),cross_den=(a*b).sum(-1),
                within_a=wa.sum(-1),within_b=wb.sum(-1),
                within_a_den=(a*(a-1)).sum(-1),within_b_den=(b*(b-1)).sum(-1))

def aggregate(raw,codes):
    return np.stack([np.sum(raw*(codes==c)[...,None],axis=-2) for c in range(4)],axis=-2)
