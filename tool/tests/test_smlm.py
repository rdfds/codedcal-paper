"""Unit tests for the count-balanced SMLM estimator (codedcal.smlm)."""
import itertools

import numpy as np

from codedcal import smlm


def test_count_balance_cancellation_identity():
    """Port of the exact enumeration check from balanced_estimator.exact():
    averaging the code-projected cubic contrast over ALL 4^3 color
    assignments of three sources equals the direct same-source triple sums,
    for both the phase-pair (within) side and the count-swapped (cross)
    side, at machine precision. This is the count-balance cancellation
    identity: both sides carry the same source-count product."""
    rng = np.random.default_rng(2609145100)
    raw = rng.normal(size=(3, 4)) + 1j * rng.normal(size=(3, 4))
    n = rng.uniform(.1, 5, size=(3, 4))
    got = []
    for colors in itertools.product(range(4), repeat=3):
        zc = np.zeros((1, 4, 4), complex)
        nc = np.zeros((1, 4, 4))
        for i, c in enumerate(colors):
            zc[0, c] += raw[i]
            nc[0, c] += n[i]
        got.append(smlm.cubic(zc, nc)[0])
    c = 0.
    w = 0.
    for i, j, k in smlm.TRIPLES:
        w += np.sum(raw[:, i] * raw[:, j].conj() * n[:, k]).real / 4
        c += (np.sum(raw[:, i] * n[:, j] * raw[:, k].conj())
              + np.sum(n[:, i] * raw[:, j] * raw[:, k].conj())).real / 8
    err = float(np.max(np.abs(np.mean(got, 0) - [c, w])))
    assert err < 1e-12


def test_window_split_components_match_pooled_cubic():
    """The window-split components, pooled arithmetically over the two
    windows, must reproduce the pooled cubic statistic exactly
    (split_components.check_matches_cubic)."""
    rng = np.random.default_rng(7)
    z = rng.normal(size=(64, 4, 4)) + 1j * rng.normal(size=(64, 4, 4))
    n = rng.uniform(0.5, 3.0, size=(64, 4, 4))
    s = smlm.components(z, n)
    assert np.allclose(s[..., :2].mean(-2), smlm.cubic(z, n))


def _synthetic_events(seed, groups, sources, tau, sigma, events_per_panel):
    """Stationary-error world: per source, one persistent displacement per
    window (variance tau^2 per axis), independent event noise sigma.
    Returns event-table columns (group, code, panel, x, weight)."""
    rng = np.random.default_rng(seed)
    G, S, K = groups, sources, events_per_panel
    theta = rng.uniform(-35, 35, (G, S))
    codes = rng.integers(0, 4, size=(G, S))
    u = rng.normal(0, tau, (G, S, 2))                  # one shift per window
    u_panel = np.repeat(u, 2, axis=-1)                 # panels 0,1=A; 2,3=B
    x = (theta[..., None, None] + u_panel[..., None]
         + rng.normal(0, sigma, (G, S, 4, K)))
    w = rng.integers(1, 6, (G, S, 4, K)).astype(float)
    group = np.broadcast_to(np.arange(G)[:, None, None, None], x.shape).ravel()
    code = np.broadcast_to(codes[..., None, None], x.shape).ravel()
    panel = np.broadcast_to(np.arange(4)[None, None, :, None], x.shape).ravel()
    return group, code, panel, x.ravel(), w.ravel()


def test_budget_recovers_total_variance_on_synthetic_events():
    """End to end from an event table: with tau=3 (T = tau^2 = 9 units^2)
    and stationary event noise, the pooled count-balanced fit lands near
    the truth and the window-split D diagnostic does not reject
    stationarity. Fixed seed; tolerances reflect the estimator's known
    finite-sample noise at this size (a few units^2 at ~4k groups)."""
    T_true = 9.0
    ev = _synthetic_events(seed=20260916, groups=4096, sources=32,
                           tau=3.0, sigma=5.0, events_per_panel=3)
    z, n, gids = smlm.aggregate_events(*ev)
    assert len(gids) == 4096
    out = smlm.budget(z, n)
    assert not out['pooled']['boundary']
    assert abs(out['T_pooled'] - T_true) < 4.0
    assert abs(out['split']['total_geometric'] - T_true) < 4.0
    # Stationary world: no strong D rejection at either frequency.
    assert np.all(np.abs(out['split']['D_z']) < 4.0)
    # Ratios must be genuine attenuations (0 < ratio < 1) at both q.
    r = np.array(out['pooled']['ratio'])
    assert np.all((r > 0) & (r < 1))


def test_aggregate_events_matches_direct_sums():
    """The event aggregator reproduces sum(w exp(iqx)) and sum(w) per
    (group, code, panel) cell exactly."""
    rng = np.random.default_rng(3)
    m = 500
    group = rng.integers(0, 7, m)
    code = rng.integers(0, 4, m)
    panel = rng.integers(0, 4, m)
    x = rng.normal(0, 20, m)
    w = rng.uniform(0.5, 2.0, m)
    q = np.array([0.12])
    z, n, gids = smlm.aggregate_events(group, code, panel, x, w, q)
    g0, c0, t0 = 2, 1, 3
    sel = (group == g0) & (code == c0) & (panel == t0)
    gi = int(np.searchsorted(gids, g0))
    assert np.isclose(n[gi, c0, t0], w[sel].sum())
    assert np.allclose(z[gi, c0, t0, 0],
                       (w[sel] * np.exp(1j * q[0] * x[sel])).sum())
