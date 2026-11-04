"""Window-split cubic components: same statistic as balanced_estimator.cubic,
resolved by which window supplies the phase pair and which panel the phase.

Under stationary source error each split ratio equals exp(-q^2 T) exactly in
population, for arbitrary heterogeneous errors and arbitrarily correlated
counts. Their ratio D is therefore a motion-free specification test, and
their geometric mean is exact under any source-uniform multiplicative change
of the error characteristic function between windows.
"""
from pathlib import Path
import sys
import numpy as np
from scipy.optimize import minimize_scalar
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'count-balanced-calibration-sep14'))
from balanced_estimator import H,Q,PERMS,old,cubic,fit_ratio,json_safe

TRIPLES_A=[(0,1,2),(0,1,3)];TRIPLES_B=[(2,3,0),(2,3,1)]

def components(z,n):
    # z,n: group,code,panel for one q. Columns per window:
    # [cross, within, orientation1, orientation2]; cross=(o1+o2)/2.
    a=np.einsum('gct,cb->gtb',z,H);b=np.einsum('gct,cb->gtb',n,H)
    def t3(x,y,w):return sum(x[:,p[0]]*y[:,p[1]]*w[:,p[2]] for p in PERMS).real/6
    out=[]
    for triples in (TRIPLES_A,TRIPLES_B):
        w=0.;o1=0.;o2=0.
        for i,j,k in triples:
            w=w+t3(a[:,i],a[:,j].conj(),b[:,k])/2
            o1=o1+t3(a[:,i],b[:,j],a[:,k].conj())/2
            o2=o2+t3(b[:,i],a[:,j],a[:,k].conj())/2
        out.append(np.stack([(o1+o2)/2,w,o1,o2],-1))
    return np.stack(out,-2) # group, window(2), component(4)

def check_matches_cubic(z,n):
    # The pooled arithmetic statistic must equal balanced_estimator.cubic.
    s=components(z,n);ref=cubic(z,n)
    got=s[...,:2].mean(-2)
    assert np.allclose(got,ref),(abs(got-ref).max())
    return True

def _fit(ratio,scale):
    t=float(minimize_scalar(lambda t:np.sum((ratio-np.exp(-scale*Q**2*t))**2),
        bounds=(-30,60),method='bounded').x)
    return t,bool(t< -29.9 or t>59.9)

def split_summary(comp):
    # comp: group,q,window,component sums over a block of groups.
    m=comp.mean(0) # q,window,component
    ra=m[:,0,0]/m[:,0,1];rb=m[:,1,0]/m[:,1,1]
    arith=(m[:,:,0].mean(-1))/(m[:,:,1].mean(-1))
    t_ar,b_ar=_fit(arith,1.)
    t_ge,b_ge=_fit(ra*rb,2.)
    # Group-level influence of log D = log ra - log rb (delta method).
    inf=((comp[:,:,0,0]-ra*comp[:,:,0,1])/(m[:,0,1]*ra)
        -(comp[:,:,1,0]-rb*comp[:,:,1,1])/(m[:,1,1]*rb))
    valid=(ra>0)&(rb>0)
    logD=np.where(valid,np.log(np.where(valid,ra/rb,1)),np.nan)
    se=inf.std(0,ddof=1)/np.sqrt(len(comp))
    orient={};orient_se={}
    for i,k in enumerate('AB'):
        r=m[:,i,2]/m[:,i,3]
        inf_o=(comp[:,:,i,2]-r*comp[:,:,i,3])/m[:,i,3]
        orient[k]=r.tolist()
        orient_se[k]=(inf_o.std(0,ddof=1)/np.sqrt(len(comp))).tolist()
    return dict(total_arithmetic=t_ar,total_geometric=t_ge,
        boundary=dict(arithmetic=b_ar,geometric=b_ge),
        ratio_A=ra.tolist(),ratio_B=rb.tolist(),
        D=(ra/rb).tolist(),logD=logD.tolist(),logD_se=se.tolist(),
        D_z=(logD/se).tolist(),
        drift_estimate_nm2=(-logD/Q**2).tolist(),
        negative_ratio=bool((~valid).any()),
        orientation_ratio=orient,orientation_se=orient_se)
