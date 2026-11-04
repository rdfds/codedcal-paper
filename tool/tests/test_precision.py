"""Precision is checked separately from exact algebraic recovery."""
import itertools
import numpy as np
import pytest
from codedcal import precision, retrofit


@pytest.fixture
def system():
    W=np.zeros((10,5),int)
    for row,bits in enumerate(itertools.combinations(range(5),3)): W[row,list(bits)]=1
    k=np.arange(10,dtype=float)+10
    R=retrofit.t_map(W,k)
    return W,k,R


def test_exact_recovery_can_fail_precision(system):
    W,k,R=system; truth=np.arange(6,dtype=float)-2.5
    C=np.array([[1.,-1,0,0,0,0]])
    V=4*(R@R.T)
    result=precision.assess(R@truth,W,k,V,max_sd=2,contrasts=C)
    assert np.allclose(result['budgets'],truth)
    assert np.allclose(result['contrast_sd'],np.sqrt(8))
    assert result['precision_rule_passed'] is False
    assert precision.assess(R@truth,W,k,V,max_sd=3,contrasts=C)['precision_rule_passed']
    assert 'blur' in result['estimand'] and 'SIMULATION ONLY' in result['status']


def test_correlated_covariance_matches_simulation(system):
    W,k,R=system; rng=np.random.default_rng(8); A=rng.normal(size=(5,5)); V=A@A.T
    result=precision.assess(np.zeros(5),W,k,V,max_sd=1000)
    draws=rng.multivariate_normal(np.zeros(5),V,size=20000)
    recovered=draws@np.linalg.pinv(R).T
    observed=np.cov(recovered,rowvar=False)
    assert np.linalg.norm(observed-result['budget_covariance'])/np.linalg.norm(observed)<.03
    assert np.max(np.abs(result['budget_covariance']@np.ones(6)))<1e-10


@pytest.mark.parametrize('bad', ['asymmetric','negative','nan','shape'])
def test_invalid_covariance_rejected(system,bad):
    W,k,R=system; V=np.eye(5)
    if bad=='asymmetric':V[0,1]=1
    if bad=='negative':V[0,0]=-1
    if bad=='nan':V[0,0]=np.nan
    if bad=='shape':V=V[:4,:4]
    with pytest.raises(ValueError):precision.assess(np.zeros(5),W,k,V,max_sd=1)


@pytest.mark.parametrize('tolerance',[0,-1,np.inf,np.nan])
def test_invalid_tolerance_rejected(system,tolerance):
    W,k,R=system
    with pytest.raises(ValueError):precision.assess(np.zeros(5),W,k,np.eye(5),max_sd=tolerance)


def test_absolute_budget_claim_rejected(system):
    W,k,R=system
    with pytest.raises(ValueError):
        precision.assess(np.zeros(5),W,k,np.eye(5),max_sd=1,contrasts=np.eye(6))


def test_singular_design_rejected():
    W=np.array([[1,1,0,0],[0,1,1,0],[0,0,1,1],[1,0,0,1]])
    with pytest.raises(ValueError):precision.assess(np.zeros(4),W,np.ones(4),np.eye(4),max_sd=1)


def test_full_rank_endpoints_can_hide_rank_loss():
    W=np.zeros((6,4),int)
    for row,bits in enumerate(itertools.combinations(range(4),2)):W[row,list(bits)]=1
    a=np.ones(6);b=np.array([20,1,1,1,1,20],float)
    assert np.linalg.matrix_rank(retrofit.t_map(W,a))==4
    assert np.linalg.matrix_rank(retrofit.t_map(W,b))==4
    result=precision.abundance_path(W,a,b)
    assert result['rank_loss_found'] and result['robustness_certified'] is False
    assert result['endpoint_inertias']==[3,2]
    assert len(result['crossings'])==1
    assert abs(result['crossings'][0]['fraction']-1/19)<1e-12
    assert result['crossings'][0]['smallest_singular']<1e-12
    reverse=precision.abundance_path(W,b,a)
    assert abs(reverse['crossings'][0]['fraction']-18/19)<1e-12


def test_no_inertia_change_does_not_certify(system):
    W,k,R=system
    result=precision.abundance_path(W,k,k*2)
    assert not result['rank_loss_found'] and not result['robustness_certified']


@pytest.mark.parametrize('bad',['weight','negative','shape','nan','uncovered'])
def test_path_invalid_inputs_rejected(system,bad):
    W,k,R=system;W=W.copy();k=k.copy()
    if bad=='weight':W[0,0]=1-W[0,0]
    if bad=='negative':k[0]=-1
    if bad=='shape':k=k[:-1]
    if bad=='nan':k[0]=np.nan
    if bad=='uncovered':W[:,0]=0
    with pytest.raises(ValueError):precision.abundance_path(W,k,k)


def test_singular_endpoint_is_retained():
    W=np.zeros((6,4),int)
    for row,bits in enumerate(itertools.combinations(range(4),2)):W[row,list(bits)]=1
    counts=np.ones(6);end=counts.copy();end[[0,5]]=2
    result=precision.abundance_path(W,counts,end)
    assert result['rank_loss_found']
    assert any(row['fraction']==1 for row in result['crossings'])
