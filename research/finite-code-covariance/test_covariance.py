import itertools
from fractions import Fraction as F
import numpy as np
import pytest
from covariance import Contrast, four_panel, kernel, gradients, cost_variance


def random_contrast(n=4):
    rng=np.random.default_rng(781)
    z=rng.normal(size=(n,4,2,2))+1j*rng.normal(size=(n,4,2,2))
    counts=rng.integers(0,5,size=(n,4))
    return four_panel(z,counts,True),z,counts


@pytest.mark.parametrize('k',[3,4,5])
def test_exact_rational_count_moments(k):
    for n in range(1,6):
        values=[]
        for codes in itertools.product(range(k),repeat=n):
            counts=[codes.count(c) for c in range(k)]
            values.append(F(k*k*sum(v**3 for v in counts)-3*k*n*sum(v*v for v in counts)+2*n**3,(k-1)*(k-2)))
        mean=sum(values)/len(values)
        var=sum((v-mean)**2 for v in values)/len(values)
        assert mean==n
        assert var==F(18*n*(n-1),k-1)+F(6*n*(n-1)*(n-2),(k-1)*(k-2))


@pytest.mark.parametrize('k',[3,4,5])
def test_exhaustive_actual_complex_covariance(k):
    c,z,n=random_contrast(4)
    out=c.evaluate(np.array(list(itertools.product(range(k),repeat=4))),k)
    analytic=c.moments(k)
    np.testing.assert_allclose(out.mean(0),analytic['mean'],rtol=2e-12,atol=2e-12)
    np.testing.assert_allclose(np.cov(out,rowvar=False,bias=True),analytic['covariance'],rtol=3e-12,atol=3e-11)
    assert np.linalg.eigvalsh(analytic['covariance']).min()>-1e-9


def test_source_contraction_against_direct_coefficients():
    c,*_=random_contrast(7)
    mu,p,q=c.direct_coefficients();m=c.moments()
    np.testing.assert_allclose(m['mean'],mu,atol=1e-12)
    np.testing.assert_allclose(m['pair_gram'],p,atol=2e-10)
    np.testing.assert_allclose(m['triple_gram'],q,atol=2e-10)


def test_shipped_hadamard_all_assignments():
    from codedcal.smlm import components,H,PERMS
    for labels in itertools.product(range(4),repeat=3):
        actual=sum(np.prod([H[labels[j],perm[j]] for j in range(3)]) for perm in PERMS)/6
        assert actual==pytest.approx(kernel(*labels,4))
    c,z,n=random_contrast(3)
    for codes in itertools.product(range(4),repeat=3):
        nz=np.zeros((1,4,4,2,2),complex);nn=np.zeros((1,4,4))
        for s,code in enumerate(codes): nz[0,code]+=z[s];nn[0,code]+=n[s]
        actual=sum(np.stack([components(nz[...,q,a],nn)[0,:,:2] for q in range(2)]) for a in range(2))/2
        got=c.evaluate(np.array(codes)).reshape(2,2,2,2)[...,0]
        np.testing.assert_allclose(actual,got,atol=3e-13)


def test_permutations_and_scaling_and_offdiagonals():
    c,z,n=random_contrast(5);m=c.moments()
    order=np.array([3,1,4,0,2]);per=four_panel(z[order],n[order],True).moments()
    np.testing.assert_allclose(per['covariance'],m['covariance'],atol=1e-10)
    scaled=four_panel(z*2,n*2,True).moments()
    np.testing.assert_allclose(scaled['covariance'],64*m['covariance'],atol=1e-8)
    assert abs(m['covariance'][0,8])>1e-3 # cross frequency is not independent
    assert abs(m['covariance'][0,4])>1e-3 # cross window is not independent
    out=c.evaluate(np.array([0,1,1,3,2]))
    relabel=c.evaluate(np.array([2,0,0,1,3]))
    np.testing.assert_allclose(out,relabel,atol=1e-12)


@pytest.mark.parametrize('n',[0,1,2])
def test_small_source_boundary(n):
    c,*_=random_contrast(n);m=c.moments()
    mu,p,q=c.direct_coefficients()
    np.testing.assert_allclose(m['mean'],mu,atol=1e-12)
    np.testing.assert_allclose(m['pair_gram'],p,atol=1e-11)
    np.testing.assert_allclose(m['triple_gram'],q,atol=1e-10)


def test_nonlinear_counterexample_and_fixed_ratio():
    ones=np.ones((3,1));c=Contrast(ones,np.array([[0,0,0]]),np.array([[1.]]))
    w=c.evaluate(np.array(list(itertools.product(range(4),repeat=3))))[:,0]
    assert set(np.round(w,10))=={-1.,3.,27.}
    assert np.var(w)==pytest.approx(42)
    d=np.linspace(.1,.9,len(w));np.testing.assert_allclose(w*d/w,d)
    # Identical C,R for every assignment: variance of C-R must be zero.
    z=np.ones((5,4,1,1),complex);n=np.ones((5,4))
    m=four_panel(z,n).moments();v=np.array([1,-1,1,-1])
    assert abs(v@m['covariance']@v)<1e-10


def test_nonuniform_assignment_breaks_diagonal_mean():
    c=Contrast(np.ones((2,1)),np.array([[0,0,0]]),np.array([[1.]]))
    # Both sources always receive code zero, instead of uniform independent codes.
    assert c.evaluate(np.array([0,0]))[0]==8
    assert c.moments()['mean'][0]==2


def test_gradients_and_resource_reversal():
    mu=np.array([2.,3.,4.,5.]);g=gradients(mu);eps=1e-6
    funcs={'arithmetic':lambda a:np.log((a[0]+a[2])/(a[1]+a[3])),
           'geometric':lambda a:.5*np.log(a[0]*a[2]/(a[1]*a[3]))}
    for name,fn in funcs.items():
        numerical=np.array([(fn(mu+eps*np.eye(4)[i])-fn(mu-eps*np.eye(4)[i]))/(2*eps) for i in range(4)])
        np.testing.assert_allclose(g[name],numerical,atol=2e-10)
    assert gradients([0,1,1,1]) is None
    assert cost_variance(.01,1,10,8)<cost_variance(.01,1,10,4)
    assert cost_variance(100,1,10,8)>cost_variance(100,1,10,4)


@pytest.mark.parametrize('bad',[True,2,3.5,-1])
def test_invalid_k(bad):
    with pytest.raises(ValueError): kernel(0,0,0,bad)


@pytest.mark.parametrize('k',[3,4,8])
def test_five_code_orbits_and_restricted_symmetrization(k):
    # One all-equal, THREE position-specific two-equal, one all-distinct.
    alpha=F(1);betas=[-alpha/F(k-1)]*3
    gamma=-(alpha+(k-1)*sum(betas))/((k-1)*(k-2))
    assert gamma==F(2,(k-1)*(k-2))
    for beta in betas:assert alpha/k+F(k-1,k)*beta==0
    rng=np.random.default_rng(845);a=rng.normal(size=(3,3,3));base=[];perturbed=[]
    # This nonsymmetric kernel is unbiased for every source coefficient tensor:
    # extra first-code term has expectation zero for each source-index pattern.
    for codes in itertools.product(range(k),repeat=3):
        t=0.;extra=0.
        for i,j,l in itertools.product(range(3),repeat=3):
            t+=a[i,j,l]*kernel(codes[i],codes[j],codes[l],k)
            extra+=a[i,j,l]*(codes[i]-(k-1)/2)
        base.append(t);perturbed.append(t+extra)
    np.testing.assert_allclose(np.mean(perturbed),np.mean(base),atol=1e-12)
    assert np.var(perturbed)>=np.var(base)-1e-12


def test_sufficient_statistic_intervals_against_shipped_api():
    from codedcal.uncertainty import fieller_set
    from codedcal.geometric_uncertainty import geometric_ratio_set
    from finite_budget import summarize_case
    rng=np.random.default_rng(172)
    for n in (4,16,64):
        for offset in (0,1,10):
            x=rng.normal(size=(n,4))+offset
            cov=np.cov(x,rowvar=False)/n
            for variant in ('arithmetic','geometric'):
                got=summarize_case(x.mean(0),cov,n,variant,.9)
                want=fieller_set(x[:,[0,2]].mean(1),x[:,[1,3]].mean(1)) if variant=='arithmetic' else geometric_ratio_set(x.reshape(n,2,2))
                assert got['kind']==want['kind']
                assert got['usable']==want['usable']
                for a,b in zip(got['intervals'],want['intervals']):
                    for key in ('lower','upper'):
                        if a[key] is None:assert b[key] is None
                        else:assert a[key]==pytest.approx(b[key],rel=1e-10,abs=1e-11)


@pytest.mark.parametrize('k',[3,4,8,16])
def test_classical_cumulant_polynomial_identity(k):
    rng=np.random.default_rng(81)
    f=rng.normal(size=(7,3));c=Contrast(f,np.array([[0,1,2]]),np.array([[1.]]))
    codes=rng.integers(k,size=7);bins=np.zeros((k,3));np.add.at(bins,codes,f)
    centered=bins-bins.mean(0)
    classical_k3=k*k*np.prod(centered,axis=1).mean()/((k-1)*(k-2))
    np.testing.assert_allclose(c.evaluate(codes,k)[0],k*classical_k3,rtol=1e-12,atol=1e-12)


def test_nonunit_common_factor_cancellation_with_unequal_counts():
    rng=np.random.default_rng(19);n=rng.integers(1,6,size=(8,4));theta=.71
    phase=np.exp(1j*np.array([0,0,theta,theta]))
    z=(n*phase)[...,None,None];c=four_panel(z,n);m=c.moments()
    d=np.cos(theta);v=np.array([1,-d,1,-d])
    assert abs(v@m['covariance']@v)<1e-9*np.linalg.norm(m['covariance'])
    codes=rng.integers(4,size=(30,8));out=c.evaluate(codes)
    np.testing.assert_allclose(out[:,[0,2]],d*out[:,[1,3]],rtol=1e-12,atol=1e-10)


def test_complex_input_rejection():
    with pytest.raises(ValueError):Contrast(np.array([[1j]]),np.array([[0,0,0]]),np.ones((1,1)))
    with pytest.raises(ValueError):four_panel(np.ones((1,4,1,1)),np.ones((1,4))+1j)


def test_exact_four_to_eight_resource_boundary():
    rng=np.random.default_rng(53)
    for a,p,q in rng.lognormal(size=(100,3)):
        assert (cost_variance(a,p,q,8)<cost_variance(a,p,q,4)) == (a<p/21+5*q/42)
    p,q=7.,11.;a=p/21+5*q/42
    assert cost_variance(a,p,q,8)==pytest.approx(cost_variance(a,p,q,4))


def test_total_covariance_and_shared_acquisition_sampling_unit():
    # Exact finite-state physical law, crossed with every code assignment.
    codes=np.array(list(itertools.product(range(4),repeat=3)))
    means=[];conditional=[];all_observations=[]
    for d in (.2,.9):
        c=Contrast(np.ones((3,1)),np.array([[0,0,0]]),np.array([[d,1.]]))
        m=c.moments();means.append(m['mean']);conditional.append(m['covariance'])
        all_observations.append(c.evaluate(codes))
    floor=np.cov(np.array(means),rowvar=False,bias=True)
    total=np.cov(np.concatenate(all_observations),rowvar=False,bias=True)
    np.testing.assert_allclose(total,floor+np.mean(conditional,axis=0),atol=1e-12)
    # More internal groups reduce code noise; they do not average a shared D.
    g=1000
    shared=floor+np.mean(conditional,axis=0)/g
    independent=(floor+np.mean(conditional,axis=0))/g
    np.testing.assert_allclose(shared-independent,(1-1/g)*floor,atol=1e-12)
    assert shared[0,0]>=floor[0,0]
    assert independent[0,0]<floor[0,0]


def test_orientation_and_phase_shipped_outputs_and_covariance():
    from codedcal.smlm import components
    from codedcal.phase_attenuation import phase_components
    _,z,n=random_contrast(3)
    c=four_panel(z,n,True,True)
    labels=np.array(list(itertools.product(range(4),repeat=3)));out=c.evaluate(labels)
    m=c.moments()
    np.testing.assert_allclose(np.cov(out,rowvar=False,bias=True),m['covariance'],atol=1e-10)
    for codes,got in zip(labels,out):
        oz=np.zeros((1,4,4,2,2),complex);on=np.zeros((1,4,4))
        for i,code in enumerate(codes):oz[0,code]+=z[i];on[0,code]+=n[i]
        real=sum(np.stack([components(oz[...,q,a],on)[0] for q in range(2)]) for a in range(2))/2
        phase=sum(np.stack([phase_components(oz[...,q,a],on)[0] for q in range(2)]) for a in range(2))/2
        expanded=got.reshape(2,2,4,2)
        np.testing.assert_allclose(expanded[...,0],real,atol=1e-12)
        np.testing.assert_allclose(expanded[:,:,[0,0,1],[0,1,0]],phase,atol=1e-12)


def test_phase_modulus_gradient_and_zero_numerator():
    from covariance import phase_log_gradient
    mu=np.array([2.,-.4,3.,.7,-1.,.5,4.,-.2]);eps=1e-6
    def target(x):
        return .5*(np.log(np.hypot(x[0],x[1]))-np.log(x[2])+np.log(np.hypot(x[4],x[5]))-np.log(x[6]))
    numerical=np.array([(target(mu+eps*np.eye(8)[i])-target(mu-eps*np.eye(8)[i]))/(2*eps) for i in range(8)])
    np.testing.assert_allclose(phase_log_gradient(mu),numerical,atol=2e-10)
    mu[:2]=0
    assert phase_log_gradient(mu) is None


def test_code_dependent_geometry_breaks_extraction_even_with_quiet_ratio():
    values=[];offsets=np.array([0.,.7,1.3,2.4])
    for codes in itertools.product(range(4),repeat=2):
        code=np.array(codes)
        z=np.repeat(np.exp(1j*offsets[code])[:,None,None,None],4,axis=1)
        out=four_panel(z,np.ones((2,4))).evaluate(code)
        values.append(out)
        np.testing.assert_allclose(out[[0,2]],out[[1,3]],atol=1e-12)
    # The source diagonal is exactly 2. Code-dependent phases leave a nonzero
    # off-diagonal expectation although the cross/reference ratio stays one.
    assert abs(np.mean(values,axis=0)[1]-2)>0.1
