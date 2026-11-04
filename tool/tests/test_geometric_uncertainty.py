"""Analytic confidence propagation and support controls for the geometric target."""
import json

import numpy as np
import pytest
from scipy.stats import t

from codedcal.geometric_uncertainty import (
    _positive_parts, _product_root, geometric_attenuation_uncertainty, geometric_ratio_set,
)
from codedcal.uncertainty import _interval


def test_dependent_windows_match_student_rectangle_and_log_transform():
    # Windows share the same fluctuations: independence between them is false.
    x = np.array([.6, .7, .8, .9, 1.])
    groups = np.stack([np.column_stack([x, np.ones(5)]),
                       np.column_stack([x/2, np.ones(5)])], axis=1)
    report = geometric_attenuation_uncertainty(groups[:, None], [.1])
    row = report['per_frequency'][0]
    radius = t.ppf(.9875, 4)*x.std(ddof=1)/np.sqrt(5)
    expected = [(x.mean()-radius)/np.sqrt(2), (x.mean()+radius)/np.sqrt(2)]
    assert row['usable'] and report['window_confidence_level'] == .975
    np.testing.assert_allclose([row['intervals'][0][k] for k in ('lower','upper')], expected)
    np.testing.assert_allclose([row['gaussian_equivalent_T_intervals'][0][k] for k in ('lower','upper')],
                               -np.log(expected[::-1])/.01)
    assert not report['physical_assumptions_verified']


def test_source_uniform_kernel_change_cancels_without_rejecting_window_above_one():
    attenuation, factor = .9, .7
    changed = np.array([[factor*attenuation, 1.], [factor*attenuation, factor**2]])
    row = geometric_ratio_set(np.broadcast_to(changed, (32,2,2)))
    assert row['ratio'] == pytest.approx(attenuation)
    assert row['windows'][1]['ratio'] > 1 and row['usable']
    arithmetic = changed[:,0].sum()/changed[:,1].sum()
    assert abs(arithmetic-attenuation) > .01


@pytest.mark.parametrize('windows,kind', [
    ([[-1.,1.],[-1.,1.]], 'empty'),
    ([[-1.,-1.],[-1.,-1.]], 'bounded'),
    ([[0.,0.],[0.,0.]], 'unbounded_positive'),
])
def test_degenerate_windows_never_earn_a_usable_positive_product(windows, kind):
    row = geometric_ratio_set(np.broadcast_to(windows, (16,2,2)))
    assert row['kind'] == kind and row['ratio'] is None and not row['usable']
    assert len(row['windows']) == 2
    json.dumps(row, allow_nan=False)


def test_weak_reference_preserves_full_disjoint_window_sets():
    y = np.array([-.9,-.9,1.,1.])
    a = np.stack([np.column_stack([np.ones(4),y])]*2,axis=1)
    row = geometric_ratio_set(a)
    assert all(w['kind'] == 'disjoint' for w in row['windows'])
    assert all(w['intervals'][0]['lower'] is None for w in row['windows'])
    assert row['ratio'] > 0 and not row['usable']
    assert row['kind'] == 'unbounded_positive' and row['intervals'][0]['upper'] is None


def test_interval_union_propagation_retains_open_zero_and_disconnected_support():
    full = [_interval(upper=1., upper_closed=True), _interval(lower=2., lower_closed=True)]
    assert _product_root(_positive_parts(full), _positive_parts(full)) == [_interval(0.,None)]
    left = [_interval(1.,1.,True,True)]
    right = [_interval(1.,4.,True,False), _interval(9.,16.,True,True)]
    assert _product_root(left,right) == [_interval(1.,2.,True,False), _interval(3.,4.,True,True)]


def test_product_positive_subset_cannot_rescue_an_unusable_window():
    a = np.array([[[1.,1.],[1.,1.]], [[-1.,1.],[1.,1.]],
                  [[1.,1.],[1.,1.]], [[1.,1.],[1.,1.]]])
    row = geometric_ratio_set(a)
    assert row['intervals'][0]['lower'] == 0. and not row['intervals'][0]['lower_closed']
    assert row['kind'] == 'bounded' and not row['usable']


def test_above_one_is_retained_as_negative_gaussian_equivalent_set():
    a = np.broadcast_to([[1.2,1.],[1.3,1.]], (8,2,2))
    row = geometric_attenuation_uncertainty(a[:,None],[.1])['per_frequency'][0]
    assert row['usable'] and row['unit_attenuation_excluded']
    assert not row['positive_motion_resolved'] and not row['attenuation_model_compatible']
    assert row['gaussian_equivalent_T_intervals'][0]['upper'] < 0


@pytest.mark.parametrize('a', [np.zeros((1,1,2,2)), np.full((8,1,2,2),np.nan)])
def test_invalid_groups_remain_distinct_from_empty_confidence_sets(a):
    row = geometric_attenuation_uncertainty(a,[.1])['per_frequency'][0]
    assert row['kind'] == 'invalid' and row['intervals'] is None
    assert row['gaussian_equivalent_T_intervals'] is None and not row['usable']
    json.dumps(row,allow_nan=False)


@pytest.mark.parametrize('confidence', [0.,1.,np.nan,np.nextafter(1.,0.)])
def test_invalid_or_unrepresentable_confidence_is_explicit(confidence):
    with pytest.raises(ValueError,match='confidence_level'):
        geometric_ratio_set(np.ones((4,2,2)),confidence)


@pytest.mark.parametrize('scale', [1e-180,1e180])
def test_extreme_common_scaling_keeps_window_support_and_product(scale):
    rng = np.random.default_rng(71)
    a = rng.uniform(.7,1.,(32,2,2))
    before, after = geometric_ratio_set(a), geometric_ratio_set(a*scale)
    assert before['kind'] == after['kind'] and before['usable'] == after['usable']
    assert before['ratio'] == pytest.approx(after['ratio'])
    np.testing.assert_allclose([before['intervals'][0][k] for k in ('lower','upper')],
                               [after['intervals'][0][k] for k in ('lower','upper')])


@pytest.mark.parametrize('q', [[],[0.],[np.nan],[[.1]]])
def test_bad_frequency_does_not_silently_change_the_target(q):
    with pytest.raises(ValueError,match='q must'):
        geometric_attenuation_uncertainty(np.ones((8,1,2,2)),q)


def test_window_shape_and_complex_contributions_are_rejected():
    with pytest.raises(ValueError,match='shape'):
        geometric_ratio_set(np.ones((8,2,3)))
    with pytest.raises(ValueError,match='real'):
        geometric_ratio_set(np.ones((8,2,2),complex))
