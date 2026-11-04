"""Analytic pre-outcome controls; no preparation performance examined."""
import inspect
import itertools
import numpy as np
from scipy.stats import t as student_t
import validation as V
from common import save, sha


def main():
    results={}; rng=np.random.default_rng(261009431)
    # Exhaust ALL 4^3 permanent assignments, arbitrary complex observations/counts.
    z=rng.normal(size=(3,4,2,4))+1j*rng.normal(size=(3,4,2,4)); n=rng.uniform(.1,10,(3,4))
    codeviews=[]
    for codes in itertools.product(range(4),repeat=3):
        oz,on=V.aggregate_observations(z,n,np.array(codes));codeviews.append(V.estimate(oz[None],on[None])[0])
    exact=V.diagonal(z,n)
    error=float(np.max(abs(np.mean(codeviews,axis=0)-exact)))
    assert np.allclose(np.mean(codeviews,axis=0),exact,atol=1e-11,rtol=1e-12)
    results['exact_64_assignment_source_extraction_max_abs']=error
    # Persistent motion and multiplicity: finite deterministic identity, not a CLT.
    counts=np.array([[3.,8.,4.,7.]])
    position=np.array([0.,0.,2.,2.])
    phase=np.exp(1j*position[None,:,None,None]*V.Q[None,None,:,None]*V.U[None,None,None,:,0])
    d=V.diagonal(phase*counts[:,:,None,None],counts)
    expected=np.mean(np.cos(2*V.Q[:,None]*V.U[None,:,0]),axis=1)
    assert np.allclose(d[:,:,0].sum(1)/d[:,:,1].sum(1),expected,atol=1e-14)
    assert np.allclose(V.diagonal(phase*counts[:,:,None,None]*3,counts*3),27*d)
    results['persistent_displacement_and_correlated_burst_identity']=True
    # Closed-form Gaussian + exact finite residual pool, stable noise with unequal counts.
    res=np.array([[-2.,0.],[2.,0.],[-2.,0.],[2.,0.],[-2.,0.],[2.,0.],[-2.,0.],[2.,0.]])
    g=dict(n=np.array([8]),source=np.zeros(8,int),df=np.array([4]),cov=np.array([[[8.,0.],[0.,.01]]]),
           res=res,lc=np.array([[3.,8.,4.,7.]]),factor=np.array([1.]))
    phi=V.kernel(g,'gaussian')
    predicted=np.exp(-.5*V.Q[None,:,None]**2*(8*V.U[None,None,:,0]**2+.01*V.U[None,None,:,1]**2))
    assert np.allclose(phi,predicted,atol=1e-14)
    residual=V.kernel(g,'residual')
    assert np.allclose(residual,np.cos(2*np.sqrt(2)*V.Q[None,:,None]*V.U[None,None,:,0]),atol=1e-14)
    for model in V.MODELS:
        for world in ('zero','homogeneous','group_common','heterogeneous','noise_shift'):
            truth,_=V.truth_group(g,model,world)
            physical=truth[0,:,:,0].sum(1)/truth[0,:,:,1].sum(1)
            assert np.allclose(physical,np.exp(-V.Q**2*(0 if world=='zero' else 3.)))
            if world!='noise_shift': assert np.allclose(truth[0],truth[1])
            else: assert np.max(abs(truth[0]-truth[1]))>0
    results['Gaussian_and_exact_residual_kernel_truth']=True
    # Fieller membership must agree with its untransformed studentized contrast.
    for j in range(100):
        x,y=rng.normal(size=(2,32));y+=j/20
        target=rng.normal();f=V.fieller(x,y,target)
        contrast=x-target*y
        expected=contrast.mean()**2<=student_t.ppf(.975,len(x)-1)**2*contrast.var(ddof=1)/len(x)
        assert f['covers']==bool(expected)
    # Degenerate finite/unbounded sets are retained, not bounded by a physical prior.
    x=np.ones(20);y=np.ones(20);f=V.fieller(x,y,1.)
    assert f['kind']=='bounded' and f['usable'] and f['covers'] and f['lower']==f['upper']==1.
    x=np.tile([-1.,1.],10);y=x.copy();f=V.fieller(x,y,1.)
    assert f['kind']=='all_real' and f['covers'] and not f['usable']
    # Source-label firewall in the estimator's public signature.
    assert list(inspect.signature(V.estimate).parameters)==['observed_z','observed_n']
    assert 'source' not in inspect.getsource(V.estimate).split('"""')[-1]
    results['Fieller_contrast_membership_100_cases_and_degenerate_sets']=True
    results['aggregate_only_estimator_interface']=True
    results['all_pass']=True
    results['validation_code_sha256']=sha(V.P/'validation.py')
    results['controls_code_sha256']=sha(__file__)
    save(V.P/'controls.json',results)
    print('All analytic controls passed',flush=True)
    return results

if __name__=='__main__': main()
