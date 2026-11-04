import numpy as np
import pytest
from codedcal.workflow_gain import anchor_relative_gains


def test_common_multiplier_cannot_change_anchored_intensity_or_threshold():
    image = np.array([[.3, .02], [.4, .01]])
    gains = np.array([2., 4.])
    reference = np.array([1., 2.])
    # The historical unanchored norm cutoff changes under an arbitrary gain gauge.
    assert not np.array_equal(np.linalg.norm(image/gains[:, None], axis=0)>.1,
                              np.linalg.norm(image/(100*gains[:, None]), axis=0)>.1)
    a, _ = anchor_relative_gains(gains, reference)
    b, _ = anchor_relative_gains(100*gains, reference)
    np.testing.assert_allclose(a, b, rtol=2e-14, atol=0)
    np.testing.assert_allclose(image/a[:, None], image/b[:, None], rtol=2e-14, atol=0)


def test_preserves_gain_ratios_and_cosine_directions():
    gains = np.array([.00001, .2, 3.])
    reference = np.array([2., 1., 5.])
    data = np.array([[1., 7.], [5., 2.], [4., 3.]])
    anchored, _ = anchor_relative_gains(gains, reference)
    before = data/gains[:, None]
    after = data/anchored[:, None]
    np.testing.assert_allclose(before/np.linalg.norm(before, axis=0), after/np.linalg.norm(after, axis=0), rtol=1e-14)
    np.testing.assert_allclose(anchored/anchored[0], gains/gains[0], rtol=1e-14)
    assert abs(np.log(anchored/reference).mean()) < 1e-14


def test_reference_sets_units_and_inputs_are_not_modified():
    gains = np.array([.3, 6.])
    reference = np.array([1., 2.])
    saved = gains.copy()
    a, report = anchor_relative_gains(gains, reference)
    b, _ = anchor_relative_gains(gains, 10*reference)
    np.testing.assert_array_equal(gains, saved)
    np.testing.assert_allclose(b, 10*a, rtol=1e-14)
    assert report['physical_gain_identified'] is False
    assert report['relative_gain_stability_established'] is False


@pytest.mark.parametrize('bad', [[], [0., 1.], [-1., 2.], [np.nan], [np.inf], [[1., 2.]]])
def test_rejects_invalid_or_incompatible_gains(bad):
    with pytest.raises(ValueError):
        anchor_relative_gains(bad, [1., 1.])
    with pytest.raises(ValueError):
        anchor_relative_gains([1., 1.], bad)


def test_extreme_common_scale_does_not_require_product_or_ratio_in_float():
    a, report = anchor_relative_gains([1e-300, 1e-300], [1e300, 1e300])
    np.testing.assert_allclose(a, [1e300, 1e300], rtol=2e-13)
    assert report['common_gain_multiplier'] is None


def test_positive_but_mismatched_reference_shape_is_rejected():
    with pytest.raises(ValueError, match='identical shapes'):
        anchor_relative_gains([1.], [1., 1.])


def test_reports_unrepresentable_gauge_instead_of_clipping():
    with pytest.raises(ValueError, match='without clipping'):
        anchor_relative_gains([1e-300, 1e300], [1e300, 1e300])


def test_new_gain_contract_cannot_reuse_an_old_isolated_workflow_cache_key():
    from types import SimpleNamespace
    from codedcal.workflow_isolated import cache_key
    case = SimpleNamespace(images=np.ones((2, 1, 2, 2)), support=np.ones((2, 1, 2, 2)),
        codebook=np.eye(2), rounds=np.array([0, 0]), colors=np.array([0, 1]),
        class_ids=np.array([0, 1]), blank_classes=np.array([], int), metadata={})
    old = dict(decoder='source_cosid', source_native_scales=[1., 2.])
    revised = dict(old, source_gain_gauge='initial_source_scaling_geometric_mean_v1')
    assert cache_key(case, old, 'barcode_model') != cache_key(case, revised, 'barcode_model')
