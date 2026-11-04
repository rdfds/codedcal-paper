"""Fiducial-free cross-round registration for coded FISH (MERFISH-class).

Assignment-free plane-pair phase correlation on unlabeled spot tables:
per-pair shifts with information weights, a robust node solve over the
plane graph, per-round trace + per-color (chromatic) offsets with a QC
gate, and comparison against the dataset's own fiducial-bead oracle.

Provenance (numerics copied verbatim; module-level grid constants became
(bins, bw) arguments, values unchanged for each source's use):
  hist_fft, pair_shift, solve_nodes, jackknife color-offset SE
      <- outputs/codedfish-port-sep16/real_tms.py
  fov_trace assembly (per-color mean-histogram subtraction, outlier-plane
  rule, round color-median trace, bead oracle, QC gate)
      <- outputs/codedfish-port-sep16/analyze_fov_scale.py (analyze_fov,
         registration part; decoding/amplitude moments not included)
  pair_shift_quadratic, bead_shift, chen2015_trace
      <- outputs/codedfish-port-sep16/real_chen2015.py
"""
import warnings

import numpy as np
from scipy.spatial import cKDTree

from . import io as cio

_GRIDS = {}


def _grid(bins, bw):
    """(kfreq, taper) for a bins x bins histogram of bin width bw.

    From real_tms.py / real_chen2015.py module constants KFREQ, QK, TAPER
    (verbatim formulas; taper sigma 1.2 px in both sources).
    """
    key = (int(bins), float(bw))
    if key not in _GRIDS:
        kfreq = np.fft.fftfreq(bins) * bins
        qk = 2 * np.pi * kfreq / (bins * bw)
        taper = np.exp(-np.add.outer(qk ** 2, qk ** 2) * 1.2 ** 2 / 2)
        _GRIDS[key] = (kfreq, taper)
    return _GRIDS[key]


def hist_fft(xy, bins=1024, bw=2.0):
    # From real_tms.py hist_fft (verbatim; range [[0, bins*bw]]^2 = [[0, 2048]]^2
    # there; real_chen2015.py uses bins=256, bw=1.0 with the same formula).
    h, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=bins,
                             range=[[0, bins * bw], [0, bins * bw]])
    return np.fft.fft2(h)


def pair_shift(Fa, Fb, bins=1024, bw=2.0, search=30., return_flags=False,
               auto_expand=False, expand_margin=1.25):
    """Sub-pixel shift d with Fa ~ content of Fb displaced by +d, plus an
    information weight (Gaussian peak amplitude / width^2).

    From real_tms.py pair_shift (verbatim; BINS/BW/KFREQ/TAPER globals ->
    (bins, bw) arguments with the same default values).

    `search` bounds the COARSE peak search, so a true shift larger than it
    cannot be found. The estimate does not then sit at the window edge —
    it lands on whatever sidelobe is tallest inside the window (a true
    40-px shift read as 24 px at search=30) — and every pair is wrong the
    same way, which keeps the node solve self-consistent and leaves the
    residual budget small. That is how a 26-px mis-registration passed QC
    with R = 0.09 px (19 Sep 2026 self-audit; outputs/wide-window-sep19).

    The detectable signature is that the GLOBAL correlation peak, which
    costs nothing extra, is taller than the best in-window peak and lies
    outside the window. return_flags=True reports it; auto_expand=True
    also refines from the global peak instead. Numerics are unchanged
    whenever the global peak is inside the window.

    `expand_margin` guards the expansion: a bare "global > in-window"
    comparison has no amplitude margin, so in a weak-peak regime (few
    shared spots) noise wins and the estimate relocks arbitrarily far
    away, turning a bounded error into an unbounded one (20 Sep 2026
    audit: on a field with 9 px of true motion, 16 of 2222 tile pairs
    relocked as far as 293 px). The global peak must now exceed the best
    in-window peak by this factor before the estimate follows it; the
    FLAG still fires on the bare comparison, so QC still sees the event.
    """
    kfreq, taper = _grid(bins, bw)
    Ppw = Fa * np.conj(Fb) * taper
    X = np.fft.ifft2(Ppw).real
    sgrid = kfreq * bw
    iw = np.flatnonzero(np.abs(sgrid) <= search)
    sub = X[np.ix_(iw, iw)]
    i, j = np.unravel_index(np.argmax(sub), sub.shape)
    d0 = np.array([sgrid[iw[i]], sgrid[iw[j]]])
    coarse = d0.copy()
    gi, gj = np.unravel_index(np.argmax(X), X.shape)
    gpk = np.array([sgrid[gi], sgrid[gj]])
    outside = bool(X[gi, gj] > sub[i, j] and np.max(np.abs(gpk)) > search)
    decisive = bool(outside and X[gi, gj] > expand_margin * max(sub[i, j], 0.))
    if decisive and auto_expand:
        d0 = gpk.copy()
    gx = d0[0] + np.arange(-6., 6.2, 0.4)
    gy = d0[1] + np.arange(-6., 6.2, 0.4)
    ex = np.exp(2j * np.pi * np.outer(kfreq, gx) / (bins * bw))
    ey = np.exp(2j * np.pi * np.outer(kfreq, gy) / (bins * bw))
    Z = np.real(ex.T @ Ppw @ ey)
    i, j = np.unravel_index(np.argmax(Z), Z.shape)
    d0 = np.array([gx[i], gy[j]])
    off = np.arange(-3., 3.15, 0.3)
    gx = d0[0] + off; gy = d0[1] + off
    ex = np.exp(2j * np.pi * np.outer(kfreq, gx) / (bins * bw))
    ey = np.exp(2j * np.pi * np.outer(kfreq, gy) / (bins * bw))
    Z = np.real(ex.T @ Ppw @ ey)
    XX, YY = np.meshgrid(off, off, indexing='ij')
    B0 = float(np.median(Z)); A0 = float(Z.max() - B0)
    i, j = np.unravel_index(np.argmax(Z), Z.shape)
    p = np.array([A0, off[i], off[j], 1.5, B0])
    for _ in range(8):
        dx = XX - p[1]; dy = YY - p[2]; r2 = dx ** 2 + dy ** 2
        e = np.exp(-r2 / (2 * p[3] ** 2))
        J = np.stack([e, p[0] * e * dx / p[3] ** 2, p[0] * e * dy / p[3] ** 2,
                      p[0] * e * r2 / p[3] ** 3, np.ones_like(e)], -1).reshape(-1, 5)
        resid = (Z - p[4] - p[0] * e).ravel()
        try:
            step = np.linalg.solve(J.T @ J + 1e-9 * np.eye(5), J.T @ resid)
        except np.linalg.LinAlgError:
            break
        p += np.clip(step, [-A0, -1, -1, -1, -A0], [A0, 1, 1, 1, A0])
        p[3] = np.clip(p[3], 0.6, 8.)
    w = max(float(p[0]) / float(p[3]) ** 2, 1e-12) if A0 > 0 else 1e-12
    d = d0 + np.clip(p[1:3], -3., 3.)
    if return_flags:
        return d, w, dict(coarse_px=[float(coarse[0]), float(coarse[1])],
                          global_peak_px=[float(gpk[0]), float(gpk[1])],
                          peak_outside_window=outside,
                          peak_outside_decisive=decisive,
                          expanded=bool(decisive and auto_expand))
    return d, w


def pair_shift_quadratic(Fa, Fb, bins=256, bw=1.0, search=6.0):
    """Sub-pixel shift via two matrix-DFT zoom stages + quadratic peak fit
    (no weight output). Used for the small dense Chen 2015 field.

    From real_chen2015.py pair_shift (verbatim; Fh globals -> (Fa, Fb)
    arguments, BINS=256/BW=1.0 defaults unchanged).
    """
    kfreq, taper = _grid(bins, bw)
    Ppw = Fa * np.conj(Fb) * taper
    X = np.fft.ifft2(Ppw).real
    sgrid = kfreq * bw
    win = np.abs(sgrid) <= search
    iw = np.flatnonzero(win)
    sub = X[np.ix_(iw, iw)]
    i, j = np.unravel_index(np.argmax(sub), sub.shape)
    d0 = np.array([sgrid[iw[i]], sgrid[iw[j]]])
    for span, step in ((1.5, 0.1), (0.21, 0.03)):
        gx = d0[0] + np.arange(-span, span + step / 2, step)
        gy = d0[1] + np.arange(-span, span + step / 2, step)
        ex = np.exp(2j * np.pi * np.outer(kfreq, gx) / (bins * bw))
        ey = np.exp(2j * np.pi * np.outer(kfreq, gy) / (bins * bw))
        Z = np.real(ex.T @ Ppw @ ey)
        if step > 0.05:
            i, j = np.unravel_index(np.argmax(Z), Z.shape)
            d0 = np.array([gx[i], gy[j]])
        else:
            XX, YY = np.meshgrid(gx - d0[0], gy - d0[1], indexing='ij')
            Xd = np.column_stack([np.ones(XX.size), XX.ravel(), YY.ravel(),
                                  XX.ravel() ** 2, YY.ravel() ** 2,
                                  (XX * YY).ravel()])
            c, *_ = np.linalg.lstsq(Xd, Z.ravel(), rcond=None)
            Hm = np.array([[2 * c[3], c[5]], [c[5], 2 * c[4]]])
            if np.linalg.det(Hm) > 1e-12 and Hm[0, 0] < 0:
                d0 = d0 + np.clip(np.linalg.solve(Hm, -c[1:3]), -0.21, 0.21)
    return d0


def solve_nodes(D, Wp, nn, pairs):
    """Robust (IRLS) weighted least-squares node positions from pairwise
    shifts D[(a,b)] ~ g_a - g_b, gauged by sum(g) = 0.

    From real_tms.py solve_nodes (verbatim; the nn=24/pairs=PAIRS defaults
    became required arguments).
    """
    w0 = np.sqrt(Wp / Wp.mean())
    A0 = np.zeros((len(pairs), nn))
    for ip, (a, b) in enumerate(pairs):
        A0[ip, a] = 1; A0[ip, b] = -1
    w = w0.copy()
    for _ in range(3):
        Aw = np.vstack([A0 * w[:, None], np.ones((1, nn))])
        rhs = np.vstack([D * w[:, None], np.zeros((1, 2))])
        sol, *_ = np.linalg.lstsq(Aw, rhs, rcond=None)
        res = np.sqrt(((D - A0 @ sol) ** 2).sum(1))
        kappa = 2.5 * max(np.median(res), 1e-9)
        w = w0 / np.sqrt(1 + (res / kappa) ** 2)
    return sol


def fov_trace(beads, spots, nr=8, nc=3, excluded_planes=((6, 1), (7, 0)),
              min_spots=200, bins=1024, bw=2.0, outlier_px=1.5,
              min_planes=16, max_outliers=4,
              subtract_mean_histogram=False, search_px=30.,
              auto_expand=False):
    """Fiducial-free trace + color offsets + QC for one FOV's reduced tables.

    ``beads``: (n, >=3) rows [round, x, y, ...] fiducials per round (used ONLY
    as the oracle for comparison; pass an empty array to skip the comparison).
    ``spots``: (m, >=4) rows [round, color_index, x, y, ...] per bit-plane.
    Coordinates in px on a [0, bins*bw)^2 field. ``excluded_planes`` are
    (round, color) planes never used for the trace (e.g. codebook-unused
    negative controls; defaults are the TMS/VZG116 unused bits).

    QC gate and outlier rule from analyze_fov_scale.py (frozen after
    calibration on TMS FOVs 001/200/400): planes deviating > ``outlier_px``
    from their round's color-median are dropped; qc_pass requires
    >= ``min_planes`` kept coding planes and <= ``max_outliers`` outliers.
    Registration numerics verbatim from analyze_fov_scale.analyze_fov;
    jackknife color-offset SEs verbatim from real_tms.analyze.
    """
    beads = np.asarray(beads, float)
    if beads.size == 0:
        beads = beads.reshape(0, 3)
    spots = np.asarray(spots, float)
    excluded = [list(p) for p in excluded_planes]
    nodes = [(r, c) for r in range(nr) for c in range(nc)]
    coding = [k for k, (r, c) in enumerate(nodes) if [r, c] not in excluded]
    # From analyze_fov_scale.py: per-plane histograms. The per-color mean-
    # histogram subtraction is OFF by default (NOSUB, the validated primary):
    # the sensitivity audit and the blind clean-room both showed it
    # manufactures ghost correlations (sensitivity.md; at scale it moves the
    # cohort median from 0.273 px to 0.135 px when removed). Pass
    # subtract_mean_histogram=True only to reproduce the historical frozen
    # variant.
    Hh = np.zeros((len(nodes), bins, bins)); nsp = []
    for k, (r, c) in enumerate(nodes):
        m = (spots[:, 0] == r) & (spots[:, 1] == c)
        xy = spots[m][:, 2:4]
        nsp.append(int(m.sum()))
        Hh[k], _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=bins,
                                     range=[[0, bins * bw], [0, bins * bw]])
    if subtract_mean_histogram:
        for c in range(nc):
            idx = [nodes.index((r, c)) for r in range(nr)]
            Hh[idx] -= Hh[idx].mean(0, keepdims=True)
    F = [np.fft.fft2(h) for h in Hh]
    keep = [k for k in coding if nsp[k] >= min_spots]
    kpairs = [(a, b) for i, a in enumerate(keep) for b in keep[i + 1:]]
    D = np.zeros((len(kpairs), 2)); Wp = np.zeros(len(kpairs))
    n_outside = n_expanded = 0
    for ip, (a, b) in enumerate(kpairs):
        D[ip], Wp[ip], fl = pair_shift(F[a], F[b], bins, bw, search_px,
                                       return_flags=True, auto_expand=auto_expand)
        n_outside += bool(fl['peak_outside_window'])
        n_expanded += bool(fl['expanded'])
    g = np.full((len(nodes), 2), np.nan)
    if keep:
        gk = solve_nodes(D, Wp, nn=len(keep),
                         pairs=[(keep.index(a), keep.index(b)) for a, b in kpairs])
        for i, k in enumerate(keep):
            g[k] = gk[i]
    # Bead oracle (analyze_fov_scale.py: coarse correlation + NN-matched
    # robust mean per round vs round 0).
    have_beads = len(beads) > 0 and (beads[:, 0] == 0).sum() > 0
    b0 = beads[beads[:, 0] == 0][:, 1:3] if have_beads else np.empty((0, 2))
    F0 = hist_fft(b0, bins, bw) if have_beads else None
    btr = [np.zeros(2)]
    for r in range(1, nr):
        br = beads[beads[:, 0] == r][:, 1:3] if have_beads else np.empty((0, 2))
        if len(br) == 0:
            btr.append(np.full(2, np.nan)); continue
        d0, _ = pair_shift(hist_fft(br, bins, bw), F0, bins, bw, search_px,
                           auto_expand=auto_expand)
        t = cKDTree(b0)
        dist, j = t.query(br - d0, distance_upper_bound=2.5)
        ok = np.isfinite(dist)
        dv = br[ok] - b0[j[ok]]
        med = np.median(dv, axis=0)
        inl = np.abs(dv - med).max(1) < 1.0
        btr.append(dv[inl].mean(0) if inl.sum() >= 3 else np.full(2, np.nan))
    btr = np.array(btr)
    if not have_beads:
        btr[0] = np.nan
    g3 = g.reshape(nr, nc, 2).copy()
    with warnings.catch_warnings():
        # Planes excluded or below min_spots are all-NaN rows; the nan-
        # aggregations below handle them by design (same behavior as the
        # source scripts) — silence only that benign RuntimeWarning.
        warnings.simplefilter('ignore', RuntimeWarning)
        # Within-round outlier planes (> outlier_px from the round color-median).
        rmed = np.nanmedian(g3, 1, keepdims=True)
        outl = np.nanmax(np.abs(g3 - rmed), -1) > outlier_px
        n_outl = int(np.nansum(outl))
        g3[outl] = np.nan
        ground = np.nanmedian(g3, 1); sig = ground - ground[0]
        diff = sig - btr
        rms = float(np.sqrt(np.nanmean(diff ** 2))) if np.isfinite(diff).any() else float('nan')
        gc = g3 - np.nanmean(g3, 1, keepdims=True)
        coloff = np.nanmean(gc, 0) - np.nanmean(gc, 0)[0]
        # Jackknife-over-rounds SE for the color offsets (real_tms.py, verbatim).
        jk = np.stack([np.nanmean(np.delete(gc, r, 0), 0)
                       - np.nanmean(np.delete(gc, r, 0), 0)[0]
                       for r in range(nr)])
        col_se = np.sqrt((nr - 1) / nr * np.nansum((jk - np.nanmean(jk, 0)) ** 2, 0))
    # QC completeness (Addendum/review fix): a PASS must mean a USABLE
    # trace — every required round estimated, counted AFTER outlier
    # removal, not merely enough planes entering the graph.
    rounds_missing = [int(r) for r in range(nr)
                      if not np.isfinite(ground[r]).all()]
    n_usable = int(np.isfinite(g3[..., 0]).sum())
    nw = float(np.median([Wp[ip] / max(nsp[a] * nsp[b], 1)
                          for ip, (a, b) in enumerate(kpairs)])) if len(Wp) else 0.
    # Consistency of the pairwise estimates with the node solution. A
    # uniformly capped search stays self-consistent, so this does NOT
    # catch that failure (the search-window flags above do); it catches
    # partial breakdowns, where some pairs lock onto a different peak.
    resid = [float(np.hypot(*(D[ip] - (g[a] - g[b]))))
             for ip, (a, b) in enumerate(kpairs)
             if np.isfinite(g[a]).all() and np.isfinite(g[b]).all()]
    max_abs_pair_shift = float(np.max(np.abs(D))) if len(D) else 0.
    return dict(n_spots_per_plane=nsp, n_planes_kept=len(keep),
                search_px=float(search_px),
                n_pairs_peak_outside_window=int(n_outside),
                n_pairs_search_expanded=int(n_expanded),
                search_expanded=bool(n_expanded),
                max_abs_pair_shift_px=max_abs_pair_shift,
                max_pair_residual_px=float(max(resid)) if resid else float('nan'),
                planes_kept=[list(nodes[k]) for k in keep],
                excluded_planes=excluded, n_outlier_planes=n_outl,
                median_pair_weight=float(np.median(Wp)) if len(Wp) else 0.,
                norm_weight=nw,
                n_planes_usable=n_usable,
                rounds_missing=rounds_missing,
                trace_complete=bool(not rounds_missing),
                # A pair whose correlation peak lies outside the search
                # window is not merely noisy: every such pair is wrong in the
                # same direction, so the node solve stays self-consistent and
                # the residual budget stays small (20 Sep 2026 audit). The
                # remedy is a wider --search, not a silent re-lock, so this
                # FAILS QC and says so.
                qc_pass=bool(n_usable >= min_planes
                             and n_outl <= max_outliers
                             and not rounds_missing
                             and n_outside == 0),
                trace_anchored=bool(np.isfinite(ground[0]).all()),
                trace_signal_px=np.round(sig, 3).tolist(),
                trace_bead_px=np.round(btr, 3).tolist(),
                trace_diff_px=np.round(diff, 3).tolist(),
                trace_rms_vs_bead_px=rms,
                color_offsets_px=np.round(coloff, 3).tolist(),
                color_offsets_se_px=np.round(col_se, 3).tolist())


# --------------------------------------------------------------------------
# Chen et al. 2015 MERFISH example: fiducial-free trace vs the dataset's own
# bead oracle (the one-command reproduction). From real_chen2015.py, verbatim.

def chen2015_trace(example_dir, rounds=16):
    """Run the fiducial-free trace vs bead oracle on the Chen 2015 example.

    ``example_dir`` must contain STORM_{r:02d}_01_c1_alist.bin (signal spots)
    and STORM_{r:02d}_01_c2_list.bin (fiducial beads) for r in 0..rounds-1.
    Returns the comparison dict (pre-registered criterion: RMS <= 0.2 px).
    """
    from pathlib import Path
    S = Path(example_dir)
    R = rounds
    pairs = [(r, rp) for r in range(R) for rp in range(r + 1, R)]
    sig, bead = [], []
    for r in range(R):
        sig.append(cio.insight3_xy(S / f'STORM_{r:02d}_01_c1_alist.bin'))
        bead.append(cio.insight3_xy(S / f'STORM_{r:02d}_01_c2_list.bin'))
    # Bleedthrough guard: drop signal spots at same-round bead positions, so
    # "fiducial-free" is not secretly bead-assisted (real_chen2015.py).
    guard_dropped = []
    for r in range(R):
        t = cKDTree(bead[r])
        d, _ = t.query(sig[r], distance_upper_bound=1.5)
        keepm = ~np.isfinite(d)
        guard_dropped.append(int((~keepm).sum()))
        sig[r] = sig[r][keepm]

    # Bead oracle trace (real_chen2015.py bead_shift + lstsq, verbatim).
    def bead_shift(r, rp):
        t = cKDTree(bead[rp])
        d, j = t.query(bead[r], distance_upper_bound=4.0)
        ok = np.isfinite(d)
        dv = bead[r][ok] - bead[rp][j[ok]]          # = g_r - g_rp
        med = np.median(dv, axis=0)
        inl = np.abs(dv - med).max(1) < 1.0
        return dv[inl].mean(0), int(inl.sum()), float(dv[inl].std(0).mean())

    brows, brhs, bres = [], [], []
    for r, rp in pairs:
        s, n, sd = bead_shift(r, rp)
        e = np.zeros(R); e[r] = 1; e[rp] = -1
        brows.append(e); brhs.append(s); bres.append((n, sd))
    brows.append(np.ones(R)); brhs.append(np.zeros(2))
    g_bead, *_ = np.linalg.lstsq(np.array(brows), np.array(brhs), rcond=None)

    # Assignment-free signal trace (256 px field, 1 px bins; verbatim grid).
    bins, bw = 256, 1.0
    Fh = [hist_fft(s, bins, bw) for s in sig]
    rows, rhs = [], []
    for r, rp in pairs:
        e = np.zeros(R); e[r] = 1; e[rp] = -1
        rows.append(e); rhs.append(pair_shift_quadratic(Fh[r], Fh[rp], bins, bw))
    rows.append(np.ones(R)); rhs.append(np.zeros(2))
    g_sig, *_ = np.linalg.lstsq(np.array(rows), np.array(rhs), rcond=None)

    diff = g_sig - g_bead
    rms = float(np.sqrt((diff ** 2).mean()))
    return dict(
        dataset='Chen 2015 MERFISH_Examples, FOV 01, 16 rounds',
        n_signal=[int(len(s)) for s in sig],
        n_bead=[int(len(b)) for b in bead],
        bleedthrough_guard_dropped=guard_dropped,
        bead_pair_inliers=dict(min=int(min(n for n, _ in bres)),
                               median=float(np.median([n for n, _ in bres])),
                               residual_sd_px=float(np.median([s for _, s in bres]))),
        trace_bead_px=np.round(g_bead, 3).tolist(),
        trace_signal_px=np.round(g_sig, 3).tolist(),
        trace_diff_px=np.round(diff, 3).tolist(),
        rms_diff_px=rms,
        max_abs_diff_px=float(np.abs(diff).max()),
        corr_x=float(np.corrcoef(g_sig[:, 0], g_bead[:, 0])[0, 1]),
        corr_y=float(np.corrcoef(g_sig[:, 1], g_bead[:, 1])[0, 1]),
        criterion='rms <= 0.2 px', passed=bool(rms <= 0.2),
        px_nm_note='STORM4 pixel ~167 nm (scale context only)')
