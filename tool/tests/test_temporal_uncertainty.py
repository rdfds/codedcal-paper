"""Independent analytic controls for confidence sets, not only fitted points."""
import json

import numpy as np
import pytest
from scipy.stats import t

from codedcal import smlm
from codedcal.uncertainty import _quadratic_set, attenuation_uncertainty, fieller_set


def test_fixed_reference_reduces_to_student_mean_interval():
    x = np.array([.7, .8, .9, 1., 1.1])
    result = fieller_set(x, np.ones(5))
    radius = t.ppf(.975, 4)*x.std(ddof=1)/np.sqrt(5)
    interval = result['intervals'][0]
    np.testing.assert_allclose([interval['lower'], interval['upper']],
                               [x.mean()-radius, x.mean()+radius])
    assert result['usable'] and result['contains_one']


@pytest.mark.parametrize('coefficients,kind,bounds', [
    ((1, -3, 2), 'bounded', [(1., 2.)]),
    ((-1, 3, -2), 'disjoint', [(None, 1.), (2., None)]),
    ((-1, 0, -1), 'all_real', [(None, None)]),
    ((1, 0, 1), 'empty', []),
    ((0, 2, -4), 'half_line', [(None, 2.)]),
    ((0, -2, 4), 'half_line', [(2., None)]),
    ((0, 0, 0), 'all_real', [(None, None)]),
    ((0, 0, 1), 'empty', []),
    ((1, -2, 1), 'bounded', [(1., 1.)]),
    ((-1, 2, -1), 'all_real', [(None, None)]),
])
def test_complete_confidence_set_topologies(coefficients, kind, bounds):
    got_kind, intervals = _quadratic_set(*coefficients)
    assert got_kind == kind
    assert [(i['lower'], i['upper']) for i in intervals] == bounds


def test_zero_mean_reference_does_not_become_a_finite_error_bar():
    result = fieller_set(np.ones(4), [-1., -1., 1., 1.])
    assert result['kind'] == 'disjoint'
    assert result['ratio'] is None
    assert not result['usable']
    assert result['intervals'][0]['lower'] is None
    assert result['intervals'][1]['upper'] is None
    json.dumps(result, allow_nan=False)


def test_zero_statistics_are_all_real_not_zero_uncertainty():
    result = fieller_set(np.zeros(8), np.zeros(8))
    assert result['kind'] == 'all_real'
    assert not result['usable']


def test_negative_ratios_retained_and_never_log_clipped():
    result = attenuation_uncertainty(np.array([[[-1., 1.]]]*8), [.1])
    row = result['per_frequency'][0]
    assert row['ratio'] == -1.
    assert row['intervals'][0]['upper'] == -1.
    assert row['gaussian_equivalent_T_intervals'] == []
    assert not row['usable']


def test_ratio_above_one_retains_negative_gaussian_equivalent_variance():
    result = attenuation_uncertainty(np.array([[[1.3, 1.]]]*8), [.1])
    row = result['per_frequency'][0]
    assert row['usable'] and row['unit_attenuation_excluded']
    assert not row['positive_motion_resolved']
    assert not row['attenuation_model_compatible']
    interval = row['gaussian_equivalent_T_intervals'][0]
    assert interval['lower'] < 0 and interval['upper'] < 0
    assert interval['lower'] == pytest.approx(-np.log(1.3)/.01)


def test_unbounded_ratio_transforms_to_unbounded_positive_part():
    x = np.ones(4)
    y = np.array([-1., -1., 1., 1.])
    row = attenuation_uncertainty(np.stack([x, y], -1)[:, None], [.1])['per_frequency'][0]
    assert row['kind'] == 'disjoint'
    assert row['gaussian_equivalent_T_intervals'][0]['lower'] is None
    assert row['gaussian_equivalent_T_intervals'][0]['upper'] is not None


@pytest.mark.parametrize('scale', [1e-180, 1e180])
def test_common_rescaling_preserves_sets_for_extreme_count_products(scale):
    x, y = np.array([.7, .9, 1.1, .8]), np.array([.9, 1.1, 1., 1.2])
    original, scaled = fieller_set(x, y), fieller_set(x*scale, y*scale)
    assert scaled['kind'] == original['kind']
    np.testing.assert_allclose(
        [scaled['intervals'][0][k] for k in ('lower', 'upper')],
        [original['intervals'][0][k] for k in ('lower', 'upper')])


@pytest.mark.parametrize('x,y,reason', [
    ([1.], [1.], 'fewer_than_two_groups'),
    ([1., np.nan], [1., 1.], 'nonfinite_group_statistics'),
])
def test_invalid_samples_are_explicit(x, y, reason):
    result = fieller_set(x, y)
    assert result['kind'] == 'invalid' and result['reason'] == reason
    json.dumps(result, allow_nan=False)


@pytest.mark.parametrize('confidence', [0., 1., float('nan')])
def test_invalid_confidence_rejected(confidence):
    with pytest.raises(ValueError, match='confidence_level'):
        fieller_set([1., 2.], [1., 2.], confidence)


def test_budget_keeps_historical_numerics_while_adding_uncertainty():
    rng = np.random.default_rng(44)
    z = rng.normal(size=(64, 4, 4, 2)) + 1j*rng.normal(size=(64, 4, 4, 2))
    n = rng.uniform(.5, 3., (64, 4, 4))
    a = np.stack([smlm.cubic(z[..., j], n) for j in range(2)], 1)
    historical_T, historical_pooled = smlm.fit_ratio(a)
    result = smlm.budget(z, n)
    assert result['T_pooled'] == historical_T
    assert result['pooled'] == historical_pooled
    assert not result['uncertainty']['physical_assumptions_verified']
