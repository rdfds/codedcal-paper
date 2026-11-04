import itertools,sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parent))
from core import (COEFF,NORM,LABELS,aggregate,coded_components,collapsed,covariance,
                  realize,population,coherence,witness_vector,collision_aware)
from inference import witness_test,hoeffding_rectangle,t_witness_diagnostic


def sample(m,seed=29):
    r=np.random.default_rng(seed);n=r.integers(1,7,(m,4)).astype(float)
    z=n*r.uniform(.1,1,(m,4))*np.exp(1j*r.uniform(-np.pi,np.pi,(m,4)))
    return z,n


def all_components(z,n):
    codes=np.array(list(itertools.product(range(4),repeat=len(z))))
    onehot=codes[:,None,:]==np.arange(4)[None,:,None]
    zz=np.einsum('gks,st->gkt',onehot,z);nn=np.einsum('gks,st->gkt',onehot,n)
    return codes,coded_components(zz,nn)


@pytest.mark.parametrize('order',[2,3,4])
def test_canonical_kernels_and_exact_norms(order):
    vals=[]
    for prefix in itertools.product(range(4),repeat=order-1):
        row=[COEFF[len(set(prefix+(c,)))] for c in range(4)]
        assert abs(sum(row))<1e-14;vals.extend(row)
    assert abs(np.mean(np.square(vals))-NORM[order])<1e-14


@pytest.mark.parametrize('sources',[1,2,3,4])
def test_conditional_covariance_matches_all_assignments(sources):
    z,n=sample(sources);codes,actual=all_components(z,n);model=covariance(z,n)
    np.testing.assert_allclose(actual.mean(0),collapsed(z,n).sum(0),rtol=2e-12,atol=2e-10)
    np.testing.assert_allclose(actual.mean(0),model['mean'],rtol=2e-12,atol=2e-10)
    np.testing.assert_allclose(np.cov(actual.T,ddof=0),model['covariance'],rtol=3e-12,atol=2e-8)
    for c,want in zip(codes,actual):np.testing.assert_allclose(realize(z,n,c),want,rtol=3e-12,atol=3e-9)
    assert np.linalg.eigvalsh(model['covariance']).min()>-1e-7


def test_original_hadamard_implementation_agrees():
    old=Path(__file__).resolve().parents[2]/'outputs/code-time-identification-sep14'
    sys.path.insert(0,str(old));from matched_code_moments import diagonal_moments
    z,n=sample(3);codes,actual=all_components(z,n)
    for i,c in enumerate(codes):
        zz,nn=aggregate(z,n,c);_,_,h4=diagonal_moments(np.concatenate([zz,zz.conj(),nn],axis=1))
        ca=(h4(0,9,6,11)+h4(0,9,10,7)+h4(8,1,6,11)+h4(8,1,10,7))/4
        want=[ca.real,ca.imag,h4(0,5,10,11).real,h4(8,9,2,7).real]
        np.testing.assert_allclose(actual[i],want,rtol=3e-12,atol=3e-10)


def test_arbitrary_changing_kernels_one_sided_and_common_factor_equality():
    r=np.random.default_rng(2910)
    for _ in range(1000):
        n=r.integers(1,8,(7,4));pa=r.uniform(.05,.4,7)*np.exp(1j*r.uniform(-3,3,7))
        pb=r.uniform(.05,1,7)*np.exp(1j*r.uniform(-3,3,7));m=r.uniform(0,1)*np.exp(1j*r.uniform(-3,3))
        assert coherence(population(n,pa,pb,m))<=abs(m)+1e-12
        c=2*np.exp(1.3j);np.testing.assert_allclose(coherence(population(n,pa,c*pa,m)),abs(m),rtol=1e-12,atol=1e-12)


def test_sharpness_mean_motion_and_heterogeneous_motion_limits():
    n=np.ones((2,4));pa=np.ones(2)
    for a in [.2,.5,.96,1.]:
        for g in [0.,a/2,a]:
            phi=np.arccos(g/a);pb=np.exp(1j*np.array([phi,-phi]))
            np.testing.assert_allclose(coherence(population(n,pa,pb,a)),g,atol=1e-14)
    assert coherence(population(n,pa,pa,np.exp(1.2j)))==pytest.approx(1.) # no mean-motion guarantee
    # Same counts, powers>=.25, unequal source motion: source-average claim fails.
    kernel=np.array([1.,.5]);mov=np.array([1.,.2]);v=np.array([(kernel**2*mov).sum(),0.,(kernel**2).sum(),(kernel**2).sum()])
    assert coherence(v)>.6 and coherence(v)<=mov.max()


def test_within_window_instability_and_correlated_error_counterexamples():
    # Per-panel kernels .25,1,1,.25: cross/reference coherence inflation25/16.
    p=np.array([.25,1.,1.,.25]);true_a=.64
    c=true_a*(p[0]+p[1])*(p[2]+p[3])/4;v=np.array([c,0.,p[0]*p[1],p[2]*p[3]])
    assert coherence(v)==pytest.approx(1.) and true_a<.96
    # An error exactly undoing physical motion leaves all observed phases fixed.
    z=np.ones((2,4),complex);n=np.ones((2,4))
    assert coherence(collapsed(z,n).sum(0))==pytest.approx(1.)
    assert np.cos(.7)<.96 # nontrivial latent physical characteristic is unobserved


def test_zero_counts_not_imputed_or_clipped():
    z,n=sample(3);z[1]=0;n[1]=0
    codes,vals=all_components(z,n)
    np.testing.assert_allclose(vals.mean(0),collapsed(z,n).sum(0),atol=2e-10)
    assert coherence([1,0,-1,1]) is None
    assert coherence([2,0,1,1])==2 # diagnostic model/sample inconsistency is retained
    with pytest.raises(ValueError):coded_components(np.ones((4,4))*2,np.ones((4,4)))


def test_collision_aware_unbiased_not_a_same_data_occupancy_filter():
    for m in [1,2,3,4]:
        z,n=sample(m);values=[collision_aware(z,n,c) for c in itertools.product(range(4),repeat=m)]
        np.testing.assert_allclose(np.mean(values,axis=0),collapsed(z,n).sum(0),atol=2e-10)
    # Source2 disappears: both assigned codes can still be distinct.
    z,n=sample(2);z[1]=0;n[1]=0
    design=collision_aware(z,n,[0,1])
    zz,nn=aggregate(z,n,[0,1]);detected=np.count_nonzero(nn.sum(1))
    assert detected==1 and np.linalg.norm(design)>0 # rejecting detected<m would bias expectation


def test_population_witness_null_and_evalue_expectation():
    n=np.ones((2,4));p=np.exp(1j*np.array([.9,-.9]));v=population(n,np.ones(2),p,.8)
    for theta in np.arange(8)*np.pi/4:
        for scale in [.5,1,2]:assert v@witness_vector(.96,theta,scale)<=1e-12
    # Conditional finite table with null coherence zero; enumerate two independent code groups.
    z=np.ones((2,4),complex);z[1,2:]=-1
    _,x=all_components(z,n)
    ev=[]
    for i,j in itertools.product(range(16),repeat=2):
        ev.append(np.exp(witness_test(x[[i,j]],np.array([16.,16.]),g=.96)['log_evalue']))
    assert np.mean(ev)<=1+1e-12


def test_joint_bounds_weak_support_and_witness_range_fail_closed():
    x=np.tile([.99,0.,1.,1.],(139,1));b=np.ones(139)
    assert not hoeffding_rectangle(x,b)['coherence_lower']>.99
    assert hoeffding_rectangle(np.zeros_like(x),b)['positive_reference_support'] is False
    with pytest.raises(ValueError):witness_test(x,b*.1)
    # Separate-reference rectangle loses precision compared with the correlated witness.
    rng=np.random.default_rng(55);r=rng.uniform(.5,1.5,139);xx=np.c_[r,np.zeros(139),r,r]
    assert t_witness_diagnostic(xx)['nominal_rejection']
    assert t_witness_diagnostic(xx)['finite_sample_guarantee'] is False


@pytest.mark.parametrize('bad',[[-1,0],[0,4],[0,1.5],[0,np.nan],[[0],[1]],['0','1']])
def test_invalid_source_codes_fail_without_dropping_sources(bad):
    z,n=sample(2)
    for function in [aggregate,realize,collision_aware]:
        with pytest.raises(ValueError):function(z,n,bad)


def test_joint_count_control_variance_is_exact_and_not_extra_data_free():
    z,n=sample(3);model=covariance(z,n,include_count=True)
    vals=np.array([realize(z,n,c,True) for c in itertools.product(range(4),repeat=3)])
    np.testing.assert_allclose(vals.mean(0),model['mean'],atol=2e-9)
    np.testing.assert_allclose(np.cov(vals.T,ddof=0),model['covariance'],rtol=1e-11,atol=2e-7)
    v=witness_vector(.96,.4,.8);w=vals[:,:4]@v;d=vals[:,4]-model['mean'][4]
    beta=(v@model['covariance'][:4,4])/model['covariance'][4,4]
    residual=w-beta*d
    expected=v@model['covariance'][:4,:4]@v-(v@model['covariance'][:4,4])**2/model['covariance'][4,4]
    assert residual.mean()==pytest.approx(w.mean(),abs=2e-9)
    assert residual.var()==pytest.approx(expected,rel=1e-11,abs=2e-7)
    assert residual.var()<=w.var()


def test_vectorized_simulation_inference_matches_independent_single_path():
    from simulation_controls import draw,batched_inference
    methods,_=draw(np.random.default_rng(20261029),3,139,.96,'stable')
    for x,capacity in methods.values():
        actual=batched_inference(x,capacity)
        for j,row in enumerate(x):
            one=witness_test(row,np.full(len(row),capacity))
            rectangle=hoeffding_rectangle(row,np.full(len(row),capacity))
            approx=t_witness_diagnostic(row)
            assert actual['mean_log_e'][j]==pytest.approx(one['log_evalue'],abs=1e-12)
            assert actual['e'][j]==one['reject_attenuation_at_most_g']
            assert actual['rectangle'][j]==(rectangle['coherence_lower']>.96)
            assert actual['positive_reference_rectangle'][j]==rectangle['positive_reference_support']
            assert actual['t_approx'][j]==approx['nominal_rejection']


def test_population_null_witness_survives_unequal_counts_and_complex_change():
    rng=np.random.default_rng(129)
    for _ in range(300):
        n=rng.integers(0,9,(5,4));pa=rng.uniform(.1,1,5)*np.exp(1j*rng.uniform(-3,3,5))
        pb=rng.uniform(.1,1,5)*np.exp(1j*rng.uniform(-3,3,5));m=.96*rng.random()*np.exp(1j*rng.uniform(-3,3))
        p=population(n,pa,pb,m)
        assert p@witness_vector(.96,rng.uniform(-3,3),rng.uniform(.2,3))<=1e-9


def test_complete_record_nonidentifiability_under_compensating_bias():
    # Both candidate physical models generate exactly the same recorded reads.
    # Changing the estimator cannot distinguish them without additional premises.
    rng=np.random.default_rng(991);sources=20;reads=64;q=.008
    x=np.broadcast_to(rng.normal(size=(sources,1,1)),(sources,4,reads)).copy()
    error=rng.normal(scale=.1,size=x.shape);bias=np.zeros_like(x)
    shifts=np.r_[np.ones(10),-np.ones(10)]*np.arccos(.8)/q
    u=np.zeros_like(x);u[:,2:]=shifts[:,None,None]
    observed=x+bias+error;alternative=x+u+bias-u+error
    np.testing.assert_allclose(observed,alternative,rtol=0,atol=2e-14)
    assert abs(np.mean(np.exp(1j*q*shifts)))==pytest.approx(.8)
    assert abs(np.mean(np.exp(1j*q*np.zeros(sources))))==1.
    # All panel sums (and therefore every coded moment) remain unchanged.
    np.testing.assert_allclose(np.exp(1j*q*observed).sum(2),
                               np.exp(1j*q*alternative).sum(2),atol=1e-13)
