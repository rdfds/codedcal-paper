"""Analytic controls for the new target-matching diagnostic."""
import numpy as np
from reference_characterize import ANGLES, Q, diagonal, phase_sums, count_weight
from simulate_reference import fast_clusters
import oracle_pipeline as O

def main():
    # A source observed at two persistent positions, with heterogeneous counts.
    counts=np.array([3,8,4,7]);xy=np.repeat([[0,0],[0,0],[2,0],[2,0]],counts,axis=0)
    pan=np.repeat(np.arange(4),counts);z,n=phase_sums(xy,pan);d=diagonal(z,n)
    expected=np.mean(np.cos(Q[:,None]*2*np.cos(ANGLES)[None,:]),axis=1)
    ratio=d[:,:,0].sum(1)/d[:,:,1].sum(1)
    assert np.allclose(ratio,expected,atol=1e-14)
    assert np.allclose(diagonal(z*np.exp(1j*Q[None,:,None]*37*np.cos(ANGLES)[None,None,:]),n),d)
    # Multiplicities change the triple weight, not a persistent-source ratio.
    z2,n2=phase_sums(np.repeat(xy,3,axis=0),np.repeat(pan,3))
    assert np.allclose(diagonal(z2,n2),27*d)
    assert count_weight(3*n)==27*count_weight(n)
    # Exact four-panel quadratic identity, including signed outputs.
    rng=np.random.default_rng(2610075);a,b,c,e=rng.normal(size=(4,100))
    C=(a+b-c-e)/2;WA=a-b;WB=c-e
    assert np.allclose((C*C-(WA*WA+WB*WB)/4)/2,
                       ((a-c)*(b-e)+(a-e)*(b-c))/4)
    # Exact partition equivalence for the sparse simulation implementation,
    # including the legacy strict distance cutoff and isolated singletons.
    fixtures=[np.zeros((0,3)),np.array([[0,0,0],[10,0,1],[19.9,0,2]]),
              *[rng.normal(size=(n,3))*20 for n in (1,5,30,100,300)]]
    for pos in fixtures:
        old=O.cluster_events(pos);new=fast_clusters(pos)
        assert np.array_equal(old[:,None]==old[None,:],new[:,None]==new[None,:])
    print('5 analytic controls and 7 strict clustering-partition controls passed')

if __name__=='__main__':main()
