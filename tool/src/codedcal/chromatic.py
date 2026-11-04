"""Assignment-free per-colour magnification, its quality gate, the
single-round (pre-decode) check, and the inter-colour displacement field.

Model. Plane k = (round r, colour c) carries a translation tau_k and its
colour's scale s_c about the field centre; for a plane pair (a, b) the shift
measured in a spatial tile centred at x is
    d_ab(x) = [tau_a - tau_b] + (s_a - s_b) (x - centre),
so the tile-to-tile variation of d_ab identifies the scale difference free
of every translation, with no decoding, assignments or fiducials. Only
cross-colour pairs carry scale; same-colour pairs are the control.

Provenance: ``estimate_scales``/``_estimate_scales``/``_tile_fft`` are
outputs/magnification-estimator-sep19/estimator.py verbatim, except that
``pair_shift`` defaults to ``codedcal.fish.pair_shift``; the 25 Sep
pilot adds optional tile/covariance diagnostics and rejects disconnected
colour graphs. ``gate`` is
validate.py ``gate`` (thresholds frozen on the development set);
``gate_single_round`` is single-round-qc-sep21 ``gate_k1``;
``displacement_field`` is zhuang-release-sep21/raw_field.py's tile fit.
Validated in registrations magnification-estimator-sep19 (Allen),
magnification-tms-sep20 (CZ Biohub TMS), barseq-crosstech-sep21 (null on
BARseq) and single-round-qc-sep21 (round 1 alone).
"""
import numpy as np

from . import fish

BINS_T = 256                     # per-tile histogram bins per axis; the tile
                                 # pixel size is set by the tile extent, so
                                 # bins * bw always equals the tile width


def _tile_fft(xy, x0, y0, step, extent, bins=BINS_T):
    """Histogram-FFT of the spots inside one tile, in tile-local coords,
    with the tile's spot CENTROID. The correlation measures the
    density-weighted mean displacement over the tile, so regressing on the
    centroid rather than the geometric centre removes a bias of a few
    percent of the slope when spots are not uniform in the tile."""
    m = ((xy[:, 0] >= x0) & (xy[:, 0] < x0 + step[0]) &
         (xy[:, 1] >= y0) & (xy[:, 1] < y0 + step[1]))
    h, _, _ = np.histogram2d(xy[m, 0] - x0, xy[m, 1] - y0, bins=bins,
                             range=[[0, extent], [0, extent]])
    cen = xy[m].mean(0) if m.sum() else np.array([x0 + step[0] / 2., y0 + step[1] / 2.])
    return np.fft.fft2(h), int(m.sum()), cen


def estimate_scales(spots, planes, colours, field_px=2048., ntile=3,
                    min_spots_tile=150, pair_shift=None, search=30.,
                    same_colour_control=True, min_ntile=2,
                    within_round_only=False, diagnostics=False, bins_per_tile=None):
    """Adaptive wrapper: sparse fields cannot fill a fine tiling (a field
    with ~400 spots per plane has no 3x3 tile above the minimum), so fall
    back to coarser tilings until the fit is determined. Coarser tiling
    means a shorter lever arm, hence larger SE, which the caller can gate on.
    bins_per_tile optionally sets histogram resolution (default 256). It
    changes numerical approximation, not physical tile support; higher values
    are experimental and do not establish better accuracy or calibrated SE.
    """
    if bins_per_tile is None:
        bins_per_tile = BINS_T
    if isinstance(bins_per_tile, (bool, np.bool_)) or not isinstance(bins_per_tile, (int, np.integer)) or bins_per_tile < 16:
        raise ValueError('bins_per_tile must be an integer >= 16')
    if pair_shift is None:
        pair_shift = fish.pair_shift
    ladder = [(nt, min_spots_tile) for nt in range(ntile, min_ntile - 1, -1)]
    ladder += [(min_ntile, max(int(min_spots_tile / f), 25)) for f in (2.5, 5.)]
    for nt, mst in ladder:
        out = _estimate_scales(spots, planes, colours, field_px, nt,
                               mst, pair_shift, search, same_colour_control,
                               within_round_only, diagnostics, bins_per_tile)
        out['ntile'] = nt; out['min_spots_tile'] = mst
        if out.get('scales') is not None:
            return out
    return out


def _estimate_scales(spots, planes, colours, field_px=2048., ntile=3,
                     min_spots_tile=150, pair_shift=None, search=30.,
                     same_colour_control=True, within_round_only=False, diagnostics=False, bins_per_tile=None):
    """Per-colour scale relative to colours[0].

    spots: (n, >=4) [round, colour_index, x, y]; planes: {k: (round, colour)};
    colours: list of colour labels in colour-index order. `pair_shift` is
    codedcal.fish.pair_shift (injected so this file stays dependency-free).
    Returns dict with scales and diagnostics. With diagnostics=True, also
    returns shared tile measurements and legacy fit covariance in colours[0]
    gauge. This covariance is not a validated confidence interval.
    """
    if bins_per_tile is None:
        bins_per_tile = BINS_T
    nc = len(colours)
    cen = field_px / 2.
    ffts, counts, cents = {}, {}, {}
    xy_of = {}
    for k, (r, c) in planes.items():
        m = (spots[:, 0] == r) & (spots[:, 1] == c)
        xy_of[k] = spots[m][:, 2:4]
    # Tile the OCCUPIED region, not the nominal frame: tissue often covers a
    # band of the field, and a fixed grid then leaves too few populated tiles
    # (found on Allen 1199651009/520, whose spots span y 514-1532).
    allxy = np.vstack([v for v in xy_of.values() if len(v)])
    lo = np.percentile(allxy, 0.5, axis=0); hi = np.percentile(allxy, 99.5, axis=0)
    step = (hi - lo) / ntile
    extent = float(max(step))            # square histogram window per tile
    bw_t = extent / bins_per_tile
    origins = [(lo[0] + i * step[0], lo[1] + j * step[1])
               for i in range(ntile) for j in range(ntile)]
    for k, xy in xy_of.items():
        for t, (x0, y0) in enumerate(origins):
            F, n, c0 = _tile_fft(xy, x0, y0, step, extent, bins=bins_per_tile)
            ffts[(k, t)] = F; counts[(k, t)] = n; cents[(k, t)] = c0 - cen

    tile_records = []
    rows = []                      # (colour_a, colour_b, slope, weight, n_tiles)
    keys = sorted(planes)
    for ia, a in enumerate(keys):
        for b in keys[ia + 1:]:
            ca, cb = planes[a][1], planes[b][1]
            if ca == cb and not same_colour_control:
                continue
            # The chromatic term lives entirely in cross-colour pairs of the
            # SAME round, which carry no inter-round motion; restricting to
            # them makes the estimate immune to rigid inter-round jumps and to
            # the search window such jumps would otherwise demand.
            if within_round_only and ca != cb and planes[a][0] != planes[b][0]:
                continue     # same-colour pairs are kept, as the control
            d, wts, ctr = [], [], []
            for t in range(ntile * ntile):
                ctr.append((cents[(a, t)] + cents[(b, t)]) / 2.)
                if counts[(a, t)] < min_spots_tile or counts[(b, t)] < min_spots_tile:
                    d.append([np.nan, np.nan]); wts.append(0.); continue
                dd, w = pair_shift(ffts[(a, t)], ffts[(b, t)], bins_per_tile, bw_t, search)
                d.append(dd); wts.append(w)
            d = np.array(d); wts = np.array(wts); centres = np.array(ctr)
            ok = np.isfinite(d).all(1) & (wts > 0)
            if diagnostics and ca != cb:
                tile_records.append(dict(pair=[int(a), int(b)], colours=[int(ca), int(cb)],
                    tile_ids=np.flatnonzero(ok).tolist(), centres=centres[ok].tolist(),
                    shifts=d[ok].tolist(), weights=wts[ok].tolist()))
            if ok.sum() < max(3, ntile):      # need >= 3 points for a slope
                continue
            # slope of shift vs tile centre, per axis, weighted
            sl = []
            for ax in range(2):
                x = centres[ok, ax]; y = d[ok, ax]; w = wts[ok]
                xm = np.average(x, weights=w); ym = np.average(y, weights=w)
                var = np.average((x - xm) ** 2, weights=w)
                if var <= 0:
                    sl.append(np.nan); continue
                sl.append(float(np.average((x - xm) * (y - ym), weights=w) / var))
            if not np.isfinite(sl).all():
                continue
            rows.append((ca, cb, float(np.mean(sl)), float(np.sum(wts[ok])), int(ok.sum())))

    # slope(a,b) estimates (s_b - s_a); solve least squares for s (s_0 = 0 offset)
    A, y, w = [], [], []
    ctrl = []
    for ca, cb, sl, wt, nt in rows:
        if ca == cb:
            ctrl.append(sl); continue
        # pair_shift(Fa, Fb) returns d with x_a = x_b + d, so
        # d(x) = (tau_a - tau_b) + (s_ca - s_cb)(x - centre): the slope
        # estimates s_ca - s_cb.
        row = np.zeros(nc); row[ca] += 1.; row[cb] -= 1.
        A.append(row); y.append(sl); w.append(wt)
    out = dict(n_pair_rows=len(rows), n_informative=len(A),
               same_colour_control_median=float(np.median(ctrl)) if ctrl else None,
               same_colour_control_n=len(ctrl))
    if diagnostics:
        out['tile_records'] = tile_records
        out['bins_per_tile'] = int(bins_per_tile)
        out['tile_bin_width_px'] = bw_t
    if len(A) < nc:
        out['scales'] = None; return out
    A = np.array(A); y = np.array(y); w = np.array(w)
    A = np.hstack([A, np.zeros((len(A), 1))])       # free constant, absorbs gauge
    A[:, -1] = 0.
    # fix s_0 = 0 by dropping its column
    Ar = A[:, 1:nc]
    # A pseudoinverse otherwise invents scales in unobserved colour components.
    if np.linalg.matrix_rank(Ar) < nc - 1:
        out.update(scales=None, failure_reason='disconnected colour graph')
        return out
    W = np.diag(w / w.sum())
    beta, *_ = np.linalg.lstsq(Ar.T @ W @ Ar, Ar.T @ W @ y, rcond=None)
    resid = y - Ar @ beta
    dof = max(len(y) - (nc - 1), 1)
    cov = np.linalg.pinv(Ar.T @ W @ Ar) * float((resid ** 2 * (w / w.sum())).sum()) / dof
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    if diagnostics:
        full_cov = np.zeros((nc, nc))
        full_cov[1:, 1:] = cov
        out['scale_covariance'] = full_cov.tolist()
        out['covariance_colour_order'] = [str(c) for c in colours]
        out['covariance_gauge'] = str(colours[0])
        out['covariance_status'] = 'legacy residual fit; not calibrated uncertainty'
    out['scales'] = {str(colours[0]): 1.0}
    out['scale_se'] = {str(colours[0]): 0.0}
    out['scale_se_gauge'] = str(colours[0])
    out['uncertainty_status'] = 'uncalibrated residual-fit diagnostic; not a confidence interval'
    out['gate_interpretation'] = 'frozen empirical eligibility rule; not an accuracy guarantee'
    for i in range(1, nc):
        out['scales'][str(colours[i])] = float(1.0 + beta[i - 1])
        out['scale_se'][str(colours[i])] = float(se[i - 1])
    out['residual_rms_permille'] = float(np.sqrt(np.mean(resid ** 2)) * 1e3)
    out['slopes_by_colour_pair'] = {
        f'{colours[ca]}-{colours[cb]}': float(np.median([r[2] for r in rows
                                                         if r[0] == ca and r[1] == cb]))
        for ca in range(nc) for cb in range(nc)
        if any(r[0] == ca and r[1] == cb for r in rows)}
    return out


# ------------------------------------------------------------------ gate --

def gate(est, se_max=6e-4, control_max=5e-4, resid_max_permille=8.):
    """Frozen quality gate (magnification-estimator-sep19 validate.gate):
    frozen empirical eligibility rule. Residual-fit SE is uncalibrated;
    neither its magnitude nor acceptance guarantees optical accuracy."""
    if est.get('scales') is None:
        return False, 'no fit'
    se = [v for k, v in est['scale_se'].items() if v > 0]
    if not se:
        return False, 'degenerate (zero SE)'
    if max(se) > se_max:
        return False, f'scale SE {max(se):.1e} > {se_max:.0e}'
    c = est.get('same_colour_control_median')
    if c is None:
        return False, 'no same-colour control available'
    if abs(c) > control_max:
        return False, f'same-colour control {c:+.1e}'
    if est.get('residual_rms_permille', 99) > resid_max_permille:
        return False, f"residual {est['residual_rms_permille']:.1f} permille > {resid_max_permille:g}"
    return True, 'pass'


def gate_single_round(est, se_max=6e-4, resid_max_permille=8.):
    """The gate for a ONE-round estimate (single-round-qc-sep21 gate_k1): the
    same-colour control needs two rounds, so it is dropped. Weaker than
    ``gate``; with one degree of freedom the SE rests on one closure residual."""
    if est.get('scales') is None:
        return False, 'no fit'
    se = [v for k, v in est['scale_se'].items() if v > 0]
    if not se:
        return False, 'degenerate (zero SE)'
    if max(se) > se_max:
        return False, f'scale SE {max(se):.1e} > {se_max:.0e}'
    if est.get('residual_rms_permille', 99) > resid_max_permille:
        return False, f"residual {est['residual_rms_permille']:.1f} permille > {resid_max_permille:g}"
    return True, 'pass'


# ------------------------------------------------------------ utilities --

def relative_to(scales, ref):
    """Scales re-expressed relative to colour ``ref``."""
    return {k: v / scales[ref] for k, v in scales.items()}


def rescale(spots, colours, scales, centre=1023.5):
    """Undo per-colour magnification: x -> centre + (x - centre) / s_colour.
    spots: (n, >=4) [round, colour_index, x, y, ...]; scales keyed by colour
    label (as returned, or re-referenced with ``relative_to``)."""
    out = np.asarray(spots, float).copy()
    for ci, c in enumerate(colours):
        s = scales[str(c)]
        m = out[:, 1] == ci
        out[m, 2:4] = centre + (out[m, 2:4] - centre) / s
    return out


def round1_qc(spots, planes, colours, first_round=0, field_px=2048., **kw):
    """Magnification from the FIRST imaging round alone, before any barcode
    exists (single-round-qc-sep21): within-round cross-colour pairs of that
    round only, judged by ``gate_single_round``. For live QC, detect spots on
    each image ALONE; a detector that subtracts a median over rounds leaks
    later rounds into the estimate."""
    pk = {k: v for k, v in planes.items() if v[0] == first_round}
    sp = np.asarray(spots)
    sp = sp[sp[:, 0] == first_round]
    est = estimate_scales(sp, pk, [str(c) for c in colours], field_px=field_px,
                          within_round_only=True, **kw)
    ok, why = gate_single_round(est)
    est['gate_pass'], est['gate_reason'] = bool(ok), why
    return est


def displacement_field(xy_a, xy_b, field_px=2048., ntile=4, search=30.,
                       min_spots_tile=150, centre=None):
    """Inter-colour displacement field between two spot sets imaging the
    same molecules (e.g. two colours of one round): tile-wise shifts fitted
    per axis as d(x) = tau + m (x - centre). Returns tau (px at the centre),
    m (relative magnification of a vs b) and the zero-displacement point
    centre - tau / m. d is a's position minus b's."""
    c0 = field_px / 2. - 0.5 if centre is None else centre
    step = np.array([field_px / ntile] * 2); ext = field_px / ntile
    bw = ext / BINS_T
    d, w, cen = [], [], []
    for i in range(ntile):
        for j in range(ntile):
            x0, y0 = i * step[0], j * step[1]
            Fa, na, ca = _tile_fft(np.asarray(xy_a, float), x0, y0, step, ext)
            Fb, nb, cb = _tile_fft(np.asarray(xy_b, float), x0, y0, step, ext)
            if na < min_spots_tile or nb < min_spots_tile:
                continue
            dd, ww = fish.pair_shift(Fa, Fb, BINS_T, bw, search)
            if np.isfinite(dd).all() and ww > 0:
                d.append(dd); w.append(ww); cen.append((ca + cb) / 2 - c0)
    if len(d) < 6:
        return None
    d, w, cen = np.array(d), np.array(w), np.array(cen)
    tau, m = np.zeros(2), np.zeros(2)
    for ax in range(2):
        A = np.column_stack([np.ones(len(d)), cen[:, ax]])
        W = np.sqrt(w)[:, None]
        beta = np.linalg.lstsq(A * W, d[:, ax] * W[:, 0], rcond=None)[0]
        tau[ax], m[ax] = beta
    mm = float(np.mean(m))
    zero = (c0 - tau / mm).tolist() if mm != 0 else None
    return dict(tau=tau.tolist(), m=m.tolist(), n_tiles=int(len(d)), zero_point=zero)
