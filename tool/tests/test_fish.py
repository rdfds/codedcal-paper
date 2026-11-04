"""Unit tests for the coded-FISH registration pipeline (codedcal.fish)."""
import numpy as np

from codedcal import fish


def test_pair_shift_recovers_known_shift():
    """A known sub-bin shift between two point sets sharing 4000 molecules
    is recovered to well under a histogram bin (bin width 2 px)."""
    rng = np.random.default_rng(11)
    bins, bw = 512, 2.0
    ext = bins * bw
    base = rng.uniform(60, ext - 60, size=(4000, 2))
    shift = np.array([5.3, -7.9])
    a = base + rng.normal(0, 0.3, base.shape)
    b = base + shift + rng.normal(0, 0.3, base.shape)
    # Convention (real_tms.py bead oracle): pair_shift(F_b, F_a) ~ +shift
    # when b = a + shift.
    d, w = fish.pair_shift(fish.hist_fft(b, bins, bw),
                           fish.hist_fft(a, bins, bw), bins, bw)
    assert w > 0
    assert np.all(np.abs(d - shift) < 0.1)


def test_pair_shift_quadratic_recovers_known_shift():
    """The Chen-2015 variant (256 px field, quadratic peak refinement)
    recovers a known sub-pixel shift."""
    rng = np.random.default_rng(12)
    bins, bw = 256, 1.0
    base = rng.uniform(20, bins - 20, size=(400, 2))
    shift = np.array([1.37, -0.62])
    a = base + rng.normal(0, 0.15, base.shape)
    b = base + shift + rng.normal(0, 0.15, base.shape)
    d = fish.pair_shift_quadratic(fish.hist_fft(b, bins, bw),
                                  fish.hist_fft(a, bins, bw), bins, bw)
    assert np.all(np.abs(d - shift) < 0.05)


def test_solve_nodes_recovers_node_positions_exactly():
    """With consistent pairwise differences the robust node solve returns
    the zero-sum node positions exactly."""
    rng = np.random.default_rng(13)
    nn = 6
    g = rng.normal(0, 3, (nn, 2))
    g -= g.mean(0)
    pairs = [(a, b) for a in range(nn) for b in range(a + 1, nn)]
    D = np.array([g[a] - g[b] for a, b in pairs])
    Wp = rng.uniform(0.5, 2.0, len(pairs))
    sol = fish.solve_nodes(D, Wp, nn, pairs)
    assert np.allclose(sol, g, atol=1e-8)


def test_fov_trace_recovers_known_trace_and_color_offsets():
    """End to end on a synthetic FOV: known per-round drift plus a known
    per-color offset, shared molecules across planes, clutter spots, and a
    separate bead set. The fiducial-free trace must match both the truth
    and the bead oracle, and the color offset must be recovered."""
    rng = np.random.default_rng(5)
    nr, nc = 4, 2
    bins, bw = 512, 2.0
    ext = bins * bw
    # Drifts chosen so every ghost lag introduced by the per-color mean-
    # histogram subtraction is >= 14 px from every true pair lag (the real
    # datasets this pipeline was validated on drift by ~2-20 px between
    # rounds; near-coincident round positions are outside its regime).
    drift = np.array([[0., 0.], [-1.1, -14.8], [-13.7, 8.7], [12.6, 6.4]])
    coloff = np.array([[0., 0.], [0.7, -0.4]])
    base = rng.uniform(80, ext - 80, size=(2500, 2))
    rows = []
    for r in range(nr):
        for c in range(nc):
            xy = base + drift[r] + coloff[c] + rng.normal(0, 0.4, base.shape)
            clutter = rng.uniform(0, ext, size=(400, 2))
            allxy = np.vstack([xy, clutter])
            rows.append(np.column_stack([np.full(len(allxy), r),
                                         np.full(len(allxy), c), allxy]))
    spots = np.vstack(rows)
    beadbase = rng.uniform(80, ext - 80, size=(300, 2))
    beads = np.vstack([np.column_stack([np.full(len(beadbase), r),
                                        beadbase + drift[r]
                                        + rng.normal(0, 0.05, beadbase.shape)])
                       for r in range(nr)])
    out = fish.fov_trace(beads, spots, nr=nr, nc=nc, excluded_planes=(),
                         min_spots=50, bins=bins, bw=bw,
                         min_planes=6, max_outliers=2)
    assert out['qc_pass']
    assert out['n_planes_kept'] == nr * nc
    sig = np.array(out['trace_signal_px'])
    expect = drift - drift[0]
    assert np.all(np.abs(sig - expect) < 0.15)
    assert out['trace_rms_vs_bead_px'] < 0.15
    co = np.array(out['color_offsets_px'])
    assert np.all(np.abs(co[1] - (coloff[1] - coloff[0])) < 0.1)


def test_fov_trace_without_beads_reports_nan_comparison():
    """No fiducials: the trace still runs; the oracle columns are NaN."""
    rng = np.random.default_rng(6)
    nr, nc = 3, 2
    bins, bw = 256, 2.0
    ext = bins * bw
    drift = np.array([[0., 0.], [12.8, 14.5], [-14.5, 14.8]])
    base = rng.uniform(40, ext - 40, size=(1200, 2))
    rows = []
    for r in range(nr):
        for c in range(nc):
            xy = base + drift[r] + rng.normal(0, 0.3, base.shape)
            rows.append(np.column_stack([np.full(len(xy), r),
                                         np.full(len(xy), c), xy]))
    spots = np.vstack(rows)
    out = fish.fov_trace(np.empty((0, 3)), spots, nr=nr, nc=nc,
                         excluded_planes=(), min_spots=50, bins=bins, bw=bw,
                         min_planes=5, max_outliers=1)
    sig = np.array(out['trace_signal_px'])
    assert np.all(np.abs(sig - drift) < 0.15)
    assert not np.isfinite(np.array(out['trace_bead_px'], float)[1]).any()
