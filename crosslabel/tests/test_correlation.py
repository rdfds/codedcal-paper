import numpy as np

from crosslabel.correlation import EDGES, label_translations, point_pair_g, registration_offsets


def test_random_points_have_flat_cross_correlation():
    rng = np.random.default_rng(6)
    xy = rng.uniform(0, 8000, (6000, 2)); lab = rng.integers(0, 4, len(xy))
    g = point_pair_g(xy, lab, rng)['g']
    c = (EDGES[:-1] + EDGES[1:]) / 2
    assert abs(g[(c > 10) & (c < 100)].mean() - 1) < 0.1


def test_registration_offset_is_recovered(molecule_field):
    # chance neighbours shrink a plain window mean towards zero; the plane-corrected
    # centroid should recover an imposed shift, also with uneven label coverage
    xy, lab, _, rng = molecule_field([0.2, 0.8], seed=7, side=14000.0, density=300.0, efficiency=0.7)
    loc = np.repeat(xy, 6, axis=0) + rng.normal(0, 3.0, (6 * len(xy), 2))
    llab = np.repeat(lab, 6)
    # uneven coverage: thin label 2 with a gradient across the field
    keep = (llab != 2) | (rng.random(len(llab)) < 0.4 + 0.6 * loc[:, 0] / loc[:, 0].max())
    loc, llab = loc[keep], llab[keep]
    shift = np.array([2.0, -1.5])
    loc[llab == 2] += shift
    off, se, _ = registration_offsets(loc, llab, fraction=0.5)
    t, _ = label_translations(off)
    est = t[2] - t[[0, 1, 3]].mean(0)
    assert np.all(np.abs(est - shift) < 4 * se[2, [0, 1, 3]].mean(0) + 0.3)
