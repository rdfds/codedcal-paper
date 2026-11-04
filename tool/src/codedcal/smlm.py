"""Count-balanced temporal calibration for coded SMLM event tables.

Estimates the total persistent displacement variance T (common + independent,
in squared coordinate units) from four event panels forming two windows
(panels 0/1 share one persistent displacement, panels 2/3 another), using
only per-event (group, code, panel, coordinate, weight) — no source
identities and no quality metadata. Also computes the window-split D
diagnostic: under stationary source error each window's ratio equals
exp(-q^2 T) in population, so D = ratio_A / ratio_B is a motion-free
specification test and the geometric mean is drift-robust.

Provenance (numerics copied verbatim; only import paths and an optional
``q`` argument, defaulting to the original constant, were adapted):
  H              <- outputs/code-time-identification-sep14/matched_code_moments.py
  Q              <- outputs/relative-code-calibration-sep14/experiment.py
  cubic, fit_ratio, PERMS, TRIPLES
                 <- outputs/count-balanced-calibration-sep14/balanced_estimator.py
  components, split_summary, _fit, TRIPLES_A/B
                 <- outputs/stationarity-diagnostic-sep15/split_components.py

Code identities are expectation identities, not finite-sample guarantees:
finite ratios are noisy and biased; retain negative estimates, invalid
ratios and boundary hits rather than clipping them.
"""
import itertools

import numpy as np
from scipy.optimize import minimize_scalar

from .uncertainty import attenuation_uncertainty
from .geometric_uncertainty import geometric_attenuation_uncertainty

# Physical four-code projector (matched_code_moments.py).
H = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
# Primary spatial frequencies, 1/coordinate-unit (experiment.py: nm^-1 there).
Q = np.array([.10, .15])
PERMS = list(itertools.permutations(range(3)))
# Panel triples: two phase panels from one window, count panel from the other.
TRIPLES = [(0, 1, 2), (0, 1, 3), (2, 3, 0), (2, 3, 1)]
TRIPLES_A = [(0, 1, 2), (0, 1, 3)]
TRIPLES_B = [(2, 3, 0), (2, 3, 1)]


def cubic(z, n):
    # From balanced_estimator.py (verbatim).
    # group,code,panel. No source or quality metadata input.
    a = np.einsum('gct,cb->gtb', z, H); b = np.einsum('gct,cb->gtb', n, H)
    def t3(x, y, z):
        return sum(x[:, p[0]] * y[:, p[1]] * z[:, p[2]] for p in PERMS).real / 6
    cross = 0.; within = 0.
    for i, j, k in TRIPLES:
        within = within + t3(a[:, i], a[:, j].conj(), b[:, k]) / 4
        cross = cross + (t3(a[:, i], b[:, j], a[:, k].conj())
                         + t3(b[:, i], a[:, j], a[:, k].conj())) / 8
    return np.stack([cross, within], -1)


def fit_ratio(a, q=None):
    # From balanced_estimator.py (verbatim; module constant Q -> argument q).
    q = Q if q is None else np.asarray(q, float)
    m = a.mean(0); ratio = m[:, 0] / m[:, 1]
    t = float(minimize_scalar(lambda t: np.sum((ratio - np.exp(-q ** 2 * t)) ** 2),
                              bounds=(-30, 60), method='bounded').x)
    return t, dict(ratio=ratio.tolist(), reference=m[:, 1].tolist(),
                   boundary=bool(t < -29.9 or t > 59.9))


def components(z, n):
    # From split_components.py (verbatim).
    # z,n: group,code,panel for one q. Columns per window:
    # [cross, within, orientation1, orientation2]; cross=(o1+o2)/2.
    a = np.einsum('gct,cb->gtb', z, H); b = np.einsum('gct,cb->gtb', n, H)
    def t3(x, y, w):
        return sum(x[:, p[0]] * y[:, p[1]] * w[:, p[2]] for p in PERMS).real / 6
    out = []
    for triples in (TRIPLES_A, TRIPLES_B):
        w = 0.; o1 = 0.; o2 = 0.
        for i, j, k in triples:
            w = w + t3(a[:, i], a[:, j].conj(), b[:, k]) / 2
            o1 = o1 + t3(a[:, i], b[:, j], a[:, k].conj()) / 2
            o2 = o2 + t3(b[:, i], a[:, j], a[:, k].conj()) / 2
        out.append(np.stack([(o1 + o2) / 2, w, o1, o2], -1))
    return np.stack(out, -2)  # group, window(2), component(4)


def _fit(ratio, scale, q):
    # From split_components.py (verbatim; module constant Q -> argument q).
    t = float(minimize_scalar(lambda t: np.sum((ratio - np.exp(-scale * q ** 2 * t)) ** 2),
                              bounds=(-30, 60), method='bounded').x)
    return t, bool(t < -29.9 or t > 59.9)


def split_summary(comp, q=None):
    # From split_components.py (verbatim; module constant Q -> argument q).
    # comp: group,q,window,component sums over a block of groups.
    q = Q if q is None else np.asarray(q, float)
    m = comp.mean(0)  # q,window,component
    ra = m[:, 0, 0] / m[:, 0, 1]; rb = m[:, 1, 0] / m[:, 1, 1]
    arith = (m[:, :, 0].mean(-1)) / (m[:, :, 1].mean(-1))
    t_ar, b_ar = _fit(arith, 1., q)
    t_ge, b_ge = _fit(ra * rb, 2., q)
    # Group-level influence of log D = log ra - log rb (delta method).
    inf = ((comp[:, :, 0, 0] - ra * comp[:, :, 0, 1]) / (m[:, 0, 1] * ra)
           - (comp[:, :, 1, 0] - rb * comp[:, :, 1, 1]) / (m[:, 1, 1] * rb))
    valid = (ra > 0) & (rb > 0)
    logD = np.where(valid, np.log(np.where(valid, ra / rb, 1)), np.nan)
    se = inf.std(0, ddof=1) / np.sqrt(len(comp))
    orient = {}; orient_se = {}
    for i, k in enumerate('AB'):
        r = m[:, i, 2] / m[:, i, 3]
        inf_o = (comp[:, :, i, 2] - r * comp[:, :, i, 3]) / m[:, i, 3]
        orient[k] = r.tolist()
        orient_se[k] = (inf_o.std(0, ddof=1) / np.sqrt(len(comp))).tolist()
    return dict(total_arithmetic=t_ar, total_geometric=t_ge,
                boundary=dict(arithmetic=b_ar, geometric=b_ge),
                ratio_A=ra.tolist(), ratio_B=rb.tolist(),
                D=(ra / rb).tolist(), logD=logD.tolist(), logD_se=se.tolist(),
                D_z=(logD / se).tolist(),
                drift_estimate_nm2=(-logD / q ** 2).tolist(),
                negative_ratio=bool((~valid).any()),
                orientation_ratio=orient, orientation_se=orient_se)


# --------------------------------------------------------------------------
# Event-table front end (packaging glue; the estimator numerics above are
# untouched — this only builds the aggregate phase/count arrays that the
# research code built inline from its simulators).

def aggregate_events(group, code, panel, x, weight=None, q=None):
    """Aggregate an event table into (z, n) arrays for the estimators.

    Parameters: per-event integer ``group`` (0-based ok, any ints), ``code``
    in 0..3, ``panel`` in 0..3 (panels 0/1 = window A, panels 2/3 = window B),
    coordinate ``x`` (units set the units of q; with the default q, nm), and
    optional positive ``weight`` (default 1, e.g. photons or spot quality).

    Returns (z, n, group_ids): z complex (groups, 4 codes, 4 panels, nq) with
    z = sum_events w * exp(i q x); n (groups, 4, 4) with n = sum_events w
    (no q dependence), matching the research aggregation
    (balanced_estimator.py: z/n built via old.aggregate of weighted sums).
    """
    q = Q if q is None else np.atleast_1d(np.asarray(q, float))
    group, code, panel = np.asarray(group), np.asarray(code), np.asarray(panel)
    x = np.asarray(x, float)
    weight = np.ones(x.shape) if weight is None else np.asarray(weight, float)
    arrays = dict(group=group, code=code, panel=panel, x=x, weight=weight)
    if any(value.ndim != 1 or value.shape != x.shape for value in arrays.values()):
        raise ValueError('event columns must be matching one-dimensional arrays')
    for name, value in arrays.items():
        if not np.issubdtype(value.dtype, np.number) or not np.isfinite(value).all():
            raise ValueError(f'{name} must contain finite numeric values')
        if name in ('group', 'code', 'panel') and not np.issubdtype(value.dtype, np.integer):
            if np.iscomplexobj(value) or np.any(value != np.floor(value)):
                raise ValueError(f'{name} must contain integers; fractional identities are not rounded')
    code, panel = code.astype(int), panel.astype(int)
    if np.any(weight <= 0):
        raise ValueError('event weights must be positive')
    if q.ndim != 1 or not len(q) or not np.isfinite(q).all() or np.any(q <= 0):
        raise ValueError('q must contain finite positive spatial frequencies')
    if not ((code >= 0) & (code < 4)).all():
        raise ValueError('code must be in 0..3')
    if not ((panel >= 0) & (panel < 4)).all():
        raise ValueError('panel must be in 0..3 (0/1 window A, 2/3 window B)')
    gids, gidx = np.unique(group, return_inverse=True)
    z = np.zeros((len(gids), 4, 4, len(q)), complex)
    n = np.zeros((len(gids), 4, 4))
    np.add.at(n, (gidx, code, panel), weight)
    for iq, qq in enumerate(q):
        np.add.at(z[..., iq], (gidx, code, panel), weight * np.exp(1j * qq * x))
    return z, n, gids


def budget(z, n, q=None, confidence_level=.95):
    """Count-balanced error budget from aggregated (z, n).

    Adds arithmetic Fieller and geometric Bonferroni-product sets to historical
    fits, whose numerical routines and return keys are preserved. These sets
    require independent groups and do not establish the physical model.

    Returns a dict with the pooled count-balanced fit of T (cubic + fit_ratio,
    balanced_estimator.py) and the window-split D diagnostic (split_components
    .py): per-window ratios, D with delta-method z-scores, the drift-robust
    geometric-mean T, and per-q drift estimates. T is the TOTAL persistent
    displacement variance (common + independent); separating the two requires
    the additional quartic common-motion contrast, which is out of scope here.
    """
    q = Q if q is None else np.atleast_1d(np.asarray(q, float))
    if q.ndim != 1 or not len(q) or not np.isfinite(q).all() or np.any(q <= 0):
        raise ValueError('q must contain finite positive spatial frequencies')
    if not np.isfinite(confidence_level) or not 0 < confidence_level < 1:
        raise ValueError('confidence_level must be finite and between 0 and 1')
    z, n = np.asarray(z), np.asarray(n)
    if z.ndim != 4 or z.shape[1:3] != (4, 4) or n.shape != z.shape[:3]:
        raise ValueError('z must have shape (groups, 4, 4, frequencies); n must match its first three axes')
    if len(z) < 2:
        raise ValueError('at least two independent groups are required for uncertainty')
    if not np.isfinite(z).all() or not np.isfinite(n).all() or np.any(n < 0):
        raise ValueError('z and n must be finite; n must be nonnegative')
    nq = z.shape[-1]
    if nq != len(q):
        raise ValueError(f'z has {nq} frequencies but q has {len(q)}')
    a = np.stack([cubic(z[..., iq], n) for iq in range(nq)], 1)
    uncertainty = attenuation_uncertainty(a, q, confidence_level)
    # Historical fits can be undefined for a zero reference. Keep their old
    # values/boundary flags; the confidence sets determine numerical usability.
    with np.errstate(divide='ignore', invalid='ignore'):
        t, pooled = fit_ratio(a, q)
        comp = np.stack([components(z[..., iq], n) for iq in range(nq)], 1)
        split = split_summary(comp, q)
    geometric = geometric_attenuation_uncertainty(comp[..., :2], q, confidence_level)
    return dict(n_groups=int(len(z)), q=q.tolist(),
                T_pooled=t, pooled=pooled, split=split, uncertainty=uncertainty,
                geometric_uncertainty=geometric)
