"""Tests for codedcal.chromatic: the assignment-free magnification
estimator, its gates, the round-1 check and the displacement field, on
synthetic truth (generator from magnification-estimator-sep19/
test_synthetic.py: one molecule population, each plane imaging a random
subset, per-round translations, varying density, 15% plane-specific junk)."""
import numpy as np

from codedcal import chromatic

FIELD, CEN = 2048., 1024.


def make(scales=(1.0, 0.997, 0.9986), nr=8, nc=3, n=6000, seed=0, jitter=0.0,
         share=0.45, junk=0.15):
    rng = np.random.default_rng(seed)
    pool = rng.uniform(60., FIELD - 60., (n, 2))
    rows, planes = [], {}
    for r in range(nr):
        tau = rng.normal(0, 3., 2) if r else np.zeros(2)
        for c in range(nc):
            k = r * nc + c; planes[k] = (r, c)
            take = rng.random(n) < share * rng.uniform(0.7, 1.3)
            xy = pool[take]
            obs = CEN + scales[c] * (xy - CEN) + tau + rng.normal(0, jitter, xy.shape)
            if junk:
                obs = np.vstack([obs, rng.uniform(60., FIELD - 60., (int(junk * len(obs)), 2))])
            rows.append(np.column_stack([np.full(len(obs), r), np.full(len(obs), c), obs]))
    return np.vstack(rows), planes


COL = ['750', '650', '560']


def test_recovers_known_magnification_and_passes_gate():
    true = (1.0, 0.997, 0.9986)
    spots, planes = make(scales=true, jitter=0.3)
    est = chromatic.estimate_scales(spots, planes, COL, within_round_only=True)
    for c, s in zip(COL, true):
        assert abs(est['scales'][c] - s / true[0]) < 2e-4
    ok, why = chromatic.gate(est)
    assert ok, why
    assert abs(est['same_colour_control_median']) < 5e-5


def test_null_magnification_reads_null():
    spots, planes = make(scales=(1.0, 1.0, 1.0), jitter=0.3, seed=3)
    est = chromatic.estimate_scales(spots, planes, COL, within_round_only=True)
    assert all(abs(v - 1.0) < 1e-4 for v in est['scales'].values())


def test_round1_only():
    true = (1.0, 0.997, 0.9986)
    spots, planes = make(scales=true, jitter=0.3, seed=5)
    est = chromatic.round1_qc(spots, planes, COL, first_round=0)
    assert est['gate_pass'], est['gate_reason']
    rel = chromatic.relative_to(est['scales'], '560')
    for c, s in zip(COL, true):
        assert abs(rel[c] - s / true[2]) < 5e-4
    # the full gate must refuse a one-round fit: it has no same-colour control
    assert chromatic.gate(est) == (False, 'no same-colour control available')


def test_rescale_inverts_magnification():
    true = (1.0, 0.997, 0.9986)
    rng = np.random.default_rng(1)
    xy = rng.uniform(100, 1900, (500, 2))
    rows = [np.column_stack([np.zeros(500), np.full(500, c), CEN - 0.5 + s * (xy - CEN + 0.5)])
            for c, s in enumerate(true)]
    back = chromatic.rescale(np.vstack(rows), COL, dict(zip(COL, true)))
    for c in range(3):
        assert np.allclose(back[back[:, 1] == c, 2:4], xy, atol=1e-9)


def test_displacement_field_recovers_offset_and_magnification():
    rng = np.random.default_rng(2)
    xy = rng.uniform(40, 2008, (40000, 2))
    m, tau = 1.8e-3, np.array([0.6, -0.4])
    c0 = 1023.5
    a = c0 + (1 + m) * (xy - c0) + tau + rng.normal(0, 0.2, xy.shape)
    b = xy + rng.normal(0, 0.2, xy.shape)
    f = chromatic.displacement_field(a, b)
    assert f is not None and f['n_tiles'] >= 12
    assert np.allclose(f['m'], m, atol=1.5e-4)
    assert np.allclose(f['tau'], tau, atol=0.1)
