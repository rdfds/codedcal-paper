"""Regression coverage for the QC success contract and the histogram
variant (18 Sep 2026, prompted by external review).

A QC PASS must mean a COMPLETE, usable trace: every round estimated.
The mean-histogram subtraction must stay OFF by default (NOSUB is the
validated primary) and reachable explicitly.
"""
import numpy as np
import pytest
from codedcal import fish

NR, NC, FIELD = 8, 3, 2048.0
EXCLUDED = ((6, 1), (7, 0))


def synth(missing_round=None, n=600, seed=7):
    rng = np.random.default_rng(seed)
    base = rng.uniform(50, FIELD - 50, size=(n, 2))
    shift = {r: np.array([0.4 * r, -0.3 * r]) for r in range(NR)}
    rows = []
    for r in range(NR):
        if r == missing_round:
            continue
        for c in range(NC):
            if (r, c) in EXCLUDED:
                continue
            xy = base + shift[r] + rng.normal(0, 0.05, size=(n, 2))
            rows += [[r, c, p[0], p[1], 1.0] for p in xy]
    b0 = rng.uniform(50, FIELD - 50, size=(60, 2))
    beads = [[r, *(p + shift[r])] for r in range(NR) for p in b0]
    return np.array(rows), np.array(beads)


def test_complete_data_passes_qc_with_complete_trace():
    spots, beads = synth()
    out = fish.fov_trace(beads, spots, nr=NR, nc=NC)
    assert out['qc_pass'] and out['trace_complete'] and out['trace_anchored']
    assert out['rounds_missing'] == []
    assert np.isfinite(np.asarray(out['trace_signal_px'], float)).all()


@pytest.mark.parametrize('missing', [0, 4])
def test_missing_round_fails_qc(missing):
    """Regression: a missing LATER round (4) once returned qc_pass=True
    with a NaN row in the trace; a missing reference round (0) leaves it
    unanchored. Both must now FAIL."""
    spots, beads = synth(missing_round=missing)
    out = fish.fov_trace(beads, spots, nr=NR, nc=NC)
    assert not out['qc_pass']
    assert not out['trace_complete']
    assert missing in out['rounds_missing']
    if missing == 0:
        assert not out['trace_anchored']


def test_nosub_is_the_default_and_variant_is_reachable():
    spots, beads = synth()
    default = fish.fov_trace(beads, spots, nr=NR, nc=NC)
    nosub = fish.fov_trace(beads, spots, nr=NR, nc=NC,
                           subtract_mean_histogram=False)
    sub = fish.fov_trace(beads, spots, nr=NR, nc=NC,
                         subtract_mean_histogram=True)
    assert default['trace_signal_px'] == nosub['trace_signal_px']
    assert default['trace_signal_px'] != sub['trace_signal_px']
