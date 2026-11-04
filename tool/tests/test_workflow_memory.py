"""Numerical controls for blocking; no public comparative outcomes are used."""
import numpy as np
import pytest
from codedcal.workflow_memory import projected_losses, striped_native_loss, BoundedGroupMixin
from codedcal.workflow_bardensr import BardensrAdapter
from codedcal.workflow_validation import ImageCase


@pytest.mark.parametrize('colors,words', [(2, 5), (4, 17), (4, 33), (4, 134)])
@pytest.mark.parametrize('block', [32, 64])
@pytest.mark.parametrize('batch', [1, 9])
def test_blocked_candidate_losses_equal_dense_full_field(colors, words, block, batch):
    rng = np.random.default_rng(825 + words)
    frames = 4 * colors
    shifted = rng.random((frames, 77, 81), dtype=np.float32)
    b = rng.random((frames, words), dtype=np.float32)
    b /= np.sqrt(np.square(b).sum(axis=0))
    mask = np.tile(np.arange(colors) < colors//2, 4)
    inner = (57, 61)
    offsets = [(y, x) for y in [-9, -8, 0, 8, 9] for x in [-9, -8, 0, 8, 9]]
    plus = (b[mask].T @ shifted[mask].reshape(mask.sum(), -1)).reshape(words, 77, 81)
    minus = (b[~mask].T @ shifted[~mask].reshape((~mask).sum(), -1)).reshape(words, 77, 81)
    denominator = float(np.square(shifted.astype(float)).sum())
    dense = []
    for y, x in offsets:
        maximum = (plus[:, 10-y:10-y+57, 10-x:10-x+61] +
                   minus[:, 10+y:10+y+57, 10+x:10+x+61]).max(axis=0).astype(float)
        dense.append(1 - float(np.square(maximum).sum()) / denominator)
    if np.__version__ != '1.26.4':
        with pytest.raises(ValueError, match='qualified NumPy 1.26.4'):
            projected_losses(b, shifted, mask, inner, offsets, denominator, block, batch)
        return
    blocked = projected_losses(b, shifted, mask, inner, offsets, denominator, block, batch)
    # Changes to BLAS blocking must not be hidden by a permissive tolerance.
    np.testing.assert_array_equal(blocked, dense)


def test_bad_block_and_empty_candidates_are_rejected():
    a = np.ones((4, 2), np.float32)
    shifted = np.ones((4, 24, 24), np.float32)
    mask = np.array([True, False, True, False])
    for block in [0, -1, True, 1.5, 1, 8]:
        with pytest.raises(ValueError):
            projected_losses(a, shifted, mask, (4, 4), [(0, 0)], 1, block)
    with pytest.raises(ValueError):
        projected_losses(a, shifted, mask, (4, 4), [], 1)
    for batch in [0, -1, True, 1.5]:
        with pytest.raises(ValueError):
            projected_losses(a, shifted, mask, (4, 4), [(0, 0)], 1, candidate_batch=batch)


def native_modules():
    pytest.importorskip('tensorflow')
    return pytest.importorskip('bardensr.registration.lowrankregistration')


@pytest.mark.parametrize('rows', [1, 7, 16])
def test_striped_native_operator_matches_original_with_full_boundary_support(rows):
    native = native_modules()
    rng = np.random.default_rng(93)
    images = rng.random((8, 13, 17), dtype=np.float32)
    book = rng.integers(0, 2, (8, 17)).astype(float)
    book[0] = 1
    maps = rng.integers(-9, 10, (8, 2))
    original = float(native.calc_loss(images.astype(float), book, -20-maps, zero_padding=20))
    blocked = striped_native_loss(images, book, maps, rows)
    assert abs(original-blocked) < 1e-12
    with pytest.raises(ValueError):
        striped_native_loss(images, book, maps+.25, rows)


@pytest.mark.parametrize('colors,zero', [(2, False), (4, False), (2, True), (4, True)])
def test_complete_search_maps_starts_and_projected_losses_match_dense(colors, zero):
    native_modules()
    rng = np.random.default_rng(127 + colors)
    frames = 2 * colors
    images = rng.random((frames, 1, 11, 13), dtype=np.float32)
    if zero:
        images[:] = 0
    book = rng.integers(0, 2, (frames, 17)).astype(float)
    book[0] = 1
    case = ImageCase(images, np.ones_like(images), book, np.repeat(np.arange(2), colors),
                     np.tile(np.arange(colors), 2), np.arange(17), np.array([], int))
    common = np.c_[np.zeros(frames), rng.integers(-2, 3, (frames, 2))]
    class Dense(BardensrAdapter):
        def __init__(self):
            pass
    class Blocked(BoundedGroupMixin, Dense):
        pass
    # The old loss is undefined on a zero-intensity image; do not turn that
    # unsupported control into a false equivalence assertion.
    if zero:
        with pytest.raises((ValueError, ZeroDivisionError)):
            Blocked().group_fit(case, common)
        return
    old_maps, old_info = Dense().group_fit(case, common)
    new_maps, new_info = Blocked().group_fit(case, common)
    np.testing.assert_array_equal(old_maps, new_maps)
    assert [r['start'] for r in old_info['group_families']] == [r['start'] for r in new_info['group_families']]
    np.testing.assert_array_equal([r['projected_loss'] for r in old_info['group_families']],
                                  [r['projected_loss'] for r in new_info['group_families']])
    assert abs(old_info['group_loss'] - new_info['group_loss']) < 1e-12
