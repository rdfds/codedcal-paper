"""Numeric regression pins for BOTH coded-FISH histogram variants.

18 Sep 2026, closing the ledger's queued item ("regression test pinning
both fish variants to frozen cohort outputs").

`test_qc_and_variants.py` already covers the CONTRACT: NOSUB is the
default, SUB is reachable, and the two differ.  It does not pin the
numbers, so a change that silently moves the estimator -- a different
peak refinement, a changed weighting, an altered outlier rule -- passes
it.  The research conclusions depend on the numbers, not the contract:
the NOSUB/SUB split is what the atlas and TMS cohort results turn on
(115/120 at median 0.135 px under NOSUB vs 46/122 at median 0.273 under
SUB), so both variants need pins, not just the primary.

Fixture is deterministic (fixed seed, no network, no large data).  The
pins below were generated from this fixture on 18 Sep 2026; regenerate
them ONLY with a recorded reason, since a change here means the
estimator moved.

Tolerances: trace and colour offsets are emitted rounded to 3 dp, so
they are pinned at 2e-3 -- tight enough that any real algorithmic drift
(which moves these by far more) fails, loose enough to survive a
last-digit rounding flip under a different FFT/BLAS build.  The
bead-oracle RMS is unrounded and is pinned relatively at 1e-5.
"""
import numpy as np
import pytest

from codedcal import fish

NR, NC, BINS, BW = 5, 3, 256, 2.0

# --- frozen expectations (generated 18 Sep 2026 from the fixture below) ---
PINS = {
    False: dict(          # NOSUB -- the validated primary, package default
        trace=[[0.0, 0.0], [-9.442, 13.241], [11.774, -10.558],
               [-13.113, -8.323], [7.688, 12.903]],
        color=[[0.0, 0.0], [0.602, -0.367], [-0.458, 0.504]],
        rms=0.02835207179891906,
    ),
    True: dict(           # SUB -- historical variant, retained for reproduction
        trace=[[0.0, 0.0], [-9.473, 13.211], [11.754, -10.588],
               [-13.128, -8.336], [7.666, 12.894]],
        color=[[0.0, 0.0], [0.596, -0.366], [-0.457, 0.499]],
        rms=0.03458858303154887,
    ),
}


def fixture(seed=2026):
    """Five rounds x three colours, known drift and per-colour offsets.

    Drifts are kept >= ~14 px apart pairwise so the ghost lags that the
    SUB variant introduces do not collide with true pair lags (same
    constraint as test_fish.py's end-to-end case).
    """
    rng = np.random.default_rng(seed)
    ext = BINS * BW
    drift = np.array([[0., 0.], [-9.4, 13.2], [11.8, -10.6],
                      [-13.1, -8.3], [7.7, 12.9]])
    coloff = np.array([[0., 0.], [0.61, -0.37], [-0.44, 0.52]])
    base = rng.uniform(80, ext - 80, size=(1800, 2))
    rows = []
    for r in range(NR):
        for c in range(NC):
            xy = base + drift[r] + coloff[c] + rng.normal(0, 0.35, base.shape)
            clutter = rng.uniform(0, ext, size=(300, 2))
            allxy = np.vstack([xy, clutter])
            rows.append(np.column_stack([np.full(len(allxy), r),
                                         np.full(len(allxy), c), allxy]))
    spots = np.vstack(rows)
    bb = rng.uniform(80, ext - 80, size=(250, 2))
    beads = np.vstack([
        np.column_stack([np.full(len(bb), r),
                         bb + drift[r] + rng.normal(0, 0.05, bb.shape)])
        for r in range(NR)])
    return spots, beads


def run(subtract):
    spots, beads = fixture()
    return fish.fov_trace(beads, spots, nr=NR, nc=NC, excluded_planes=(),
                          min_spots=50, bins=BINS, bw=BW, min_planes=8,
                          max_outliers=3,
                          subtract_mean_histogram=subtract)


@pytest.mark.parametrize('subtract', [False, True],
                         ids=['nosub_primary', 'sub_historical'])
def test_variant_numbers_are_pinned(subtract):
    out = run(subtract)
    pin = PINS[subtract]
    assert out['qc_pass']
    np.testing.assert_allclose(
        np.asarray(out['trace_signal_px'], float),
        np.asarray(pin['trace'], float), atol=2e-3,
        err_msg='fiducial-free trace moved; the estimator changed')
    np.testing.assert_allclose(
        np.asarray(out['color_offsets_px'], float),
        np.asarray(pin['color'], float), atol=2e-3,
        err_msg='per-colour offsets moved; the estimator changed')
    np.testing.assert_allclose(
        float(out['trace_rms_vs_bead_px']), pin['rms'], rtol=1e-5,
        err_msg='bead-oracle RMS moved; the estimator changed')


def test_nosub_beats_sub_against_the_bead_oracle():
    """Directional pin on the reason NOSUB is primary.

    On this fixture NOSUB sits closer to the bead oracle than SUB
    (0.0284 vs 0.0346 px).  The ordering, not just the values, is what
    the NOSUB adoption rests on, so it is asserted separately: a change
    that preserved both numbers' magnitudes but flipped the ordering
    would be a substantive regression.
    """
    nosub = run(False)['trace_rms_vs_bead_px']
    sub = run(True)['trace_rms_vs_bead_px']
    assert nosub < sub, (f'NOSUB ({nosub:.5f}) should beat SUB ({sub:.5f}) '
                         'against the bead oracle on this fixture')


def test_default_matches_the_nosub_pin():
    """The package default must BE the pinned primary, not merely differ
    from SUB (which test_qc_and_variants.py already checks)."""
    default = fish.fov_trace(*fixture()[::-1], nr=NR, nc=NC,
                             excluded_planes=(), min_spots=50, bins=BINS,
                             bw=BW, min_planes=8, max_outliers=3)
    np.testing.assert_allclose(
        np.asarray(default['trace_signal_px'], float),
        np.asarray(PINS[False]['trace'], float), atol=2e-3)


# --- outlier path -------------------------------------------------------
# The clean fixture above leaves max_outliers and min_planes inert (a
# mutation probe on 18 Sep moved neither the trace nor the RMS), so the
# outlier rule was uncovered.  This case corrupts one plane and pins the
# ACTUAL semantics, which are worth recording because they are not the
# obvious ones: max_outliers gates the QC VERDICT, it does not change the
# estimate, and one bad plane in fifteen degrades the trace ~4x while
# still passing QC at the default setting.

# PIN HISTORY (stated, not absorbed). The corrupting offset is 53 px,
# beyond the 30-px coarse search, so these values record pairs locked onto
# a sidelobe inside the window. On 19 Sep fov_trace briefly re-estimated
# such pairs from the global correlation peak (rms 0.10967); the 20 Sep
# audit showed that silent re-locking is unbounded in weak-peak regimes
# (worst case 226 px vs 44 px without it), so expansion is now OFF by
# default and a peak outside the window FAILS QC instead. That restores
# these original values, and the corrupted fixture now reports qc_pass
# False with 14 flagged pairs.
CORRUPT_RMS = 0.11377038932785513
CORRUPT_TRACE = [[0.0, 0.0], [-9.452, 13.232], [11.771, -10.561],
                 [-13.361, -8.078], [7.687, 12.889]]


def corrupted(max_outliers):
    spots, beads = fixture()
    m = (spots[:, 0] == 3) & (spots[:, 1] == 1)
    bad = spots.copy()
    bad[m, 2:] += np.array([40., -35.])      # bogus 40 px offset on one plane
    return fish.fov_trace(beads, bad, nr=NR, nc=NC, excluded_planes=(),
                          min_spots=50, bins=BINS, bw=BW, min_planes=8,
                          max_outliers=max_outliers)


def test_corrupted_plane_degrades_trace_but_is_pinned():
    out = corrupted(3)
    np.testing.assert_allclose(float(out['trace_rms_vs_bead_px']),
                               CORRUPT_RMS, rtol=1e-5)
    np.testing.assert_allclose(np.asarray(out['trace_signal_px'], float),
                               np.asarray(CORRUPT_TRACE, float), atol=2e-3)
    # ~4x worse than the clean fixture, but still inside the 0.25 px bar
    # the cohort analyses pre-registered.
    assert float(out['trace_rms_vs_bead_px']) > 3 * PINS[False]['rms']
    assert float(out['trace_rms_vs_bead_px']) < 0.25


@pytest.mark.parametrize('max_outliers,expect_pass', [(3, False), (1, False),
                                                      (0, False)])
def test_max_outliers_gates_qc_only(max_outliers, expect_pass):
    """Documented semantics: the threshold moves the QC verdict, never the
    estimate. If a future change makes it affect the solve, this fails and
    the change gets stated rather than absorbed. Since 20 Sep every variant
    here fails QC for a separate reason - the fixture's 53-px corruption puts
    14 pair peaks outside the 30-px search window - so the estimate, not the
    verdict, is what this test now pins."""
    out = corrupted(max_outliers)
    assert out['qc_pass'] is expect_pass
    np.testing.assert_allclose(float(out['trace_rms_vs_bead_px']),
                               CORRUPT_RMS, rtol=1e-5)


def test_search_window_failure_is_detected_and_fails_qc():
    """A shift larger than the coarse search window does NOT land at the
    window edge - it lands on a sidelobe inside it, and every pair is wrong
    the same way, so the node solve stays consistent and the residual budget
    stays small. The detectable signature is the global correlation peak
    sitting outside the window; the remedy is a wider search, so this FAILS
    QC rather than being silently re-estimated (20 Sep 2026)."""
    out = corrupted(3)
    assert out['n_pairs_peak_outside_window'] == 14
    assert not out['qc_pass']                      # detected, and refused
    assert not out['search_expanded']              # no silent re-lock
    clean = fish.fov_trace(*fixture()[::-1], nr=NR, nc=NC, excluded_planes=(),
                           min_spots=50, bins=BINS, bw=BW, min_planes=8)
    assert clean['n_pairs_peak_outside_window'] == 0 and clean['qc_pass']
    np.testing.assert_allclose(float(clean['trace_rms_vs_bead_px']),
                               PINS[False]['rms'], rtol=1e-5)


def test_opt_in_expansion_is_bounded_by_a_margin():
    """auto_expand is available but OFF by default, and it only follows a
    global peak that beats the in-window peak by expand_margin. Without the
    margin, a weak-peak regime relocks arbitrarily far away (the 20 Sep audit
    measured 226 px against a 44 px bound)."""
    rng = np.random.default_rng(3)
    b = rng.uniform(100, 1900, (400, 2)); off = np.array([40., -5.])
    Fa, Fb = fish.hist_fft(b + off, 1024, 2.), fish.hist_fft(b, 1024, 2.)
    d, _, fl = fish.pair_shift(Fa, Fb, 1024, 2., 30., return_flags=True,
                               auto_expand=True)
    np.testing.assert_allclose(d, off, atol=0.1)   # genuine jump still found
    assert fl['peak_outside_window'] and fl['expanded']
    d_off = fish.pair_shift(Fa, Fb, 1024, 2., 30., auto_expand=False)[0]
    assert np.hypot(*(d_off - off)) > 10.          # bounded, and wrong
    # a huge margin refuses every expansion, so the estimate is the capped one
    d_m = fish.pair_shift(Fa, Fb, 1024, 2., 30., auto_expand=True,
                          expand_margin=1e9)[0]
    np.testing.assert_allclose(d_m, d_off, atol=1e-9)


def test_capped_search_is_self_consistent_so_residuals_cannot_see_it():
    """Why the residual budget missed it. Three planes: two coincident, one
    displaced 53 px. With the search capped below the displacement, both
    pairs involving it lock onto the SAME wrong peak, so the node solve
    reproduces the pairwise estimates almost exactly — a small residual
    budget on a grossly mis-registered field. The window flag fires."""
    rng = np.random.default_rng(7)
    xy = rng.uniform(100., 1900., (600, 2))
    off = np.array([40., -35.])                      # 53 px
    F = [fish.hist_fft(xy, BINS, BW), fish.hist_fft(xy, BINS, BW),
         fish.hist_fft(xy + off, BINS, BW)]
    pairs = [(0, 1), (0, 2), (1, 2)]
    D, W, flags = [], [], []
    for a, b in pairs:
        d, w, fl = fish.pair_shift(F[a], F[b], BINS, BW, 30.,
                                   return_flags=True, auto_expand=False)
        D.append(d); W.append(w); flags.append(fl)
    g = fish.solve_nodes(np.array(D), np.array(W), nn=3, pairs=pairs)
    resid = max(float(np.hypot(*(np.array(D)[i] - (g[a] - g[b]))))
                for i, (a, b) in enumerate(pairs))
    assert resid < 0.5                                # self-consistent
    assert np.hypot(*(g[2] - g[0] - off)) > 10.       # yet badly wrong
    assert [f['peak_outside_window'] for f in flags] == [False, True, True]
    # with expansion the same three pairs recover the true displacement
    D2 = [fish.pair_shift(F[a], F[b], BINS, BW, 30., auto_expand=True)[0]
          for a, b in pairs]
    g2 = fish.solve_nodes(np.array(D2), np.array(W), nn=3, pairs=pairs)
    np.testing.assert_allclose(g2[2] - g2[0], off, atol=0.2)
