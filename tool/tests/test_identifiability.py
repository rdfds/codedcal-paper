"""Physical-model counterexamples, not just numerical self-consistency."""
import numpy as np
import pytest
from codedcal.identifiability import assess, relative_contrasts
from codedcal import certificate, retrofit


def test_gauge_and_constructive_witness():
    A = np.array([[1., 1.]])
    result = assess(A, np.array([[1., 1.], [1., -1.]]))
    assert result['targets'][0]['identified']
    witness = np.array(result['targets'][1]['parameter_witness'])
    np.testing.assert_allclose(A @ witness, 0, atol=1e-12)
    assert abs(np.array([1., -1.]) @ witness) > 1
    np.testing.assert_allclose(np.array(result['targets'][0]['observation_weights']) @ A, [1., 1.], atol=1e-12)


def test_one_blur_observation_does_not_separate_every_budget():
    W = certificate.rm24_shell(); kappa = np.ones(len(W)); n = W.shape[1]+1
    for H in (None, np.eye(n)[:1]):
        result = retrofit.assess_blur_model(W, kappa, blur_measurements=H)
        assert not result['all_targets_identified']
        assert all(not t['identified'] for t in result['targets'])
        for t in result['targets']:
            d = np.array(t['parameter_witness']); R = retrofit.t_map(W, kappa)
            np.testing.assert_allclose(R @ (d[:n]+d[n:]), 0, atol=1e-12)
            if H is not None: np.testing.assert_allclose(H @ d[n:], 0, atol=1e-12)
            assert t['target_change'] > .1


def test_relative_blur_observations_suffice_without_absolute_level():
    W = certificate.rm24_shell(); n = W.shape[1]+1
    H = relative_contrasts(n)
    result = retrofit.assess_blur_model(W, np.ones(len(W)), blur_measurements=H)
    assert result['all_targets_identified']
    assert result['nullity'] == 2  # separate global T and blur levels
    assert max(t['reconstruction_residual'] for t in result['targets']) < 1e-10


def test_common_blur_is_a_sufficient_assumption_not_an_inference():
    W = certificate.rm24_shell(); n = W.shape[1]+1
    result = retrofit.assess_blur_model(W, np.ones(len(W)), blur_basis=np.ones((n, 1)))
    assert result['all_targets_identified']
    assert not result['physical_blur_observations_verified']


def test_constant_weight_alone_does_not_guarantee_rank():
    W = np.array([[1,1,0,0], [0,1,1,0], [0,0,1,1], [1,0,0,1]])
    assert not retrofit.identified_map(W, np.ones(4))['all_contrasts_identified']
    with pytest.raises(ValueError, match='identifiable'):
        retrofit.solve_budgets(np.zeros(4), W, np.ones(4))


@pytest.mark.parametrize('W,kappa', [([[1,0],[0,1]],[1,1]),
    ([[1,1,0]],[1]), ([[1,1]],[-1]), ([[1,1]],[np.nan]), ([[1,2]],[1])])
def test_invalid_or_unobserved_designs_refused(W, kappa):
    with pytest.raises(ValueError): retrofit.contrast_weights(W, kappa)


def test_empty_observation_model_does_not_certify_nonzero_target():
    result=assess(np.empty((0,2)), np.eye(2))
    assert not result['all_targets_identified']
    assert result['nullity']==2
