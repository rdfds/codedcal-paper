"""Source Data 71 geometric attenuation with conservative joint Fieller sets.

Port of temporal-variants-oct17/methods.py, retaining positive-set propagation
and the requirement that BOTH original window sets be usable. The arithmetic
Fieller implementation stays unchanged for historical reproduction.
"""
import numpy as np

from .uncertainty import _interval, fieller_set


def _positive_parts(intervals):
    if intervals is None:
        return None
    result = []
    for part in intervals:
        lo, hi = part['lower'], part['upper']
        if hi is not None and hi <= 0:
            continue
        result.append(_interval(max(0., lo) if lo is not None else 0., hi,
                                part['lower_closed'] if lo is not None and lo > 0 else False,
                                part['upper_closed']))
    return result


def _product_root(left, right):
    """Image of positive interval unions under sqrt(x*y), without clipping."""
    if left is None or right is None:
        return None
    result = []
    for a in left:
        for b in right:
            lo = float(np.sqrt(a['lower'])*np.sqrt(b['lower']))
            hi = (None if a['upper'] is None or b['upper'] is None else
                  float(np.sqrt(a['upper'])*np.sqrt(b['upper'])))
            result.append(_interval(lo, hi, a['lower_closed'] and b['lower_closed'],
                                    hi is not None and a['upper_closed'] and b['upper_closed']))
    result.sort(key=lambda r: r['lower'])
    merged = []
    for part in result:
        if not merged:
            merged.append(part.copy())
            continue
        last = merged[-1]
        touches = (last['upper'] is None or part['lower'] < last['upper'] or
                   (part['lower'] == last['upper'] and
                    (last['upper_closed'] or part['lower_closed'])))
        if not touches:
            merged.append(part.copy())
            continue
        if last['upper'] is None:
            continue
        if part['upper'] is None or part['upper'] > last['upper']:
            last['upper'], last['upper_closed'] = part['upper'], part['upper_closed']
        elif part['upper'] == last['upper']:
            last['upper_closed'] |= part['upper_closed']
    return merged


def _contains(parts, value):
    return any((p['lower'] is None or value > p['lower'] or
                (value == p['lower'] and p['lower_closed'])) and
               (p['upper'] is None or value < p['upper'] or
                (value == p['upper'] and p['upper_closed'])) for p in parts or [])


def geometric_ratio_set(window_statistics, confidence_level=.95):
    """Set for sqrt((E C_A/E R_A)*(E C_B/E R_B)) at one frequency.

    Input is (independent groups, two paired windows, numerator/reference).
    Each marginal uses confidence (1+confidence_level)/2. Bonferroni does not
    require independence between windows; marginal Fieller coverage is still
    approximate. Two negative ratios never earn a positive point estimate.
    """
    if not np.isfinite(confidence_level) or not 0 < confidence_level < 1:
        raise ValueError('confidence_level must be finite and between 0 and 1')
    if np.iscomplexobj(window_statistics):
        raise ValueError('window_statistics must be real group contributions')
    a = np.asarray(window_statistics, float)
    if a.ndim != 3 or a.shape[1:] != (2, 2):
        raise ValueError('window_statistics must have shape (groups, 2, 2)')
    marginal = (1+float(confidence_level))/2
    if marginal >= 1:
        raise ValueError('confidence_level is too close to 1 for two-window allocation')
    windows = [fieller_set(a[:, w, 0], a[:, w, 1], marginal) for w in range(2)]
    # Do not reinterpret a numerical overflow as an infinite statistical bound.
    for f in windows:
        finite_set = all(p[k] is None or np.isfinite(p[k])
                         for p in f['intervals'] or [] for k in ('lower', 'upper'))
        if not finite_set:
            f.update(kind='invalid', intervals=None, usable=False,
                     contains_one=None, reason='nonfinite_confidence_calculation')
        if f['ratio'] is not None and not np.isfinite(f['ratio']):
            f['ratio'] = None
    with np.errstate(over='ignore', under='ignore'):
        parts = _product_root(*(_positive_parts(f['intervals']) for f in windows))
    finite_product = all(p[k] is None or np.isfinite(p[k])
                         for p in parts or [] for k in ('lower', 'upper'))
    if not finite_product:
        parts = None
    positive_reference = all(f['positive_reference'] for f in windows)
    available = positive_reference and all(f['ratio'] is not None and f['ratio'] > 0 for f in windows)
    ratio = float(np.exp(np.mean([np.log(f['ratio']) for f in windows]))) if available else None
    usable = bool(all(f['usable'] for f in windows) and finite_product)
    kind = ('invalid' if parts is None else 'empty' if not parts else
            'disjoint' if len(parts) > 1 else
            'bounded' if parts[0]['upper'] is not None else 'unbounded_positive')
    reason = (None if usable else 'invalid_window_confidence_set' if any(f['kind'] == 'invalid' for f in windows)
              else 'nonfinite_product_confidence_set' if not finite_product
              else 'nonpositive_window_reference' if not positive_reference
              else 'unusable_window_confidence_set')
    return dict(confidence_level=float(confidence_level), window_confidence_level=marginal,
                n_groups=len(a), ratio=ratio, usable=usable, kind=kind, intervals=parts,
                positive_reference=positive_reference, contains_one=(None if parts is None else _contains(parts, 1.)),
                reason=reason, windows=windows)


def geometric_attenuation_uncertainty(group_statistics, q, confidence_level=.95):
    """Pointwise geometric sets for (groups, frequencies, windows, C/R).

    Full window sets remain visible. The geometric set uses only strictly
    positive window ratios. A positive product cannot rescue an unusable window.
    No confidence set for the historical pooled multi-frequency fit is implied.
    """
    if np.iscomplexobj(group_statistics) or np.iscomplexobj(q):
        raise ValueError('group_statistics and q must be real')
    a, q = np.asarray(group_statistics, float), np.asarray(q, float)
    if q.ndim != 1 or not len(q) or not np.isfinite(q).all() or np.any(q <= 0):
        raise ValueError('q must contain finite positive spatial frequencies')
    if a.ndim != 4 or a.shape[1:] != (len(q), 2, 2):
        raise ValueError('group_statistics must have shape (groups, len(q), 2, 2)')
    reports = []
    for j, frequency in enumerate(q):
        f = geometric_ratio_set(a[:, j], confidence_level)
        transformed = [_interval(
            float(-np.log(p['upper'])/frequency**2) if p['upper'] is not None else None,
            float(-np.log(p['lower'])/frequency**2) if p['lower'] > 0 else None,
            p['upper_closed'] if p['upper'] is not None else False,
            p['lower_closed'] if p['lower'] > 0 else False) for p in f['intervals'] or []]
        f.update(q=float(frequency), gaussian_equivalent_T_intervals=(
                     None if f['kind'] == 'invalid' else transformed),
                 unit_attenuation_excluded=bool(f['usable'] and not f['contains_one']),
                 positive_motion_resolved=bool(f['usable'] and f['intervals'][0]['upper'] < 1),
                 attenuation_model_compatible=(None if f['kind'] == 'invalid' else any(
                     p['lower'] < 1 or (p['lower'] == 1 and p['lower_closed']) for p in f['intervals'])))
        reports.append(f)
    return dict(method='pointwise conservative geometric Fieller product set, Bonferroni rectangle',
                confidence_level=float(confidence_level), window_confidence_level=(1+float(confidence_level))/2,
                interval_endpoint_convention='null lower/upper means unbounded below/above',
                target='sqrt((E C_A/E R_A)*(E C_B/E R_B)) with positive window ratios and references',
                interpretation='geometric attenuation; -log(ratio)/q^2 is a Gaussian-equivalent variance only conditionally',
                assumptions=['independent groups; paired windows may be dependent',
                             'physical source-code exogeneity and conditionally independent temporal panels',
                             'stable source kernels, or the specified source-uniform positive kernel change under common motion',
                             'common displacement law for a common variance interpretation'],
                limitations=['approximate marginal coverage; conservative joint rectangle, not exact coverage',
                             'pointwise in frequency, not simultaneous across frequencies or variants',
                             'concentrated contributions and source-dependent or code-dependent violations can defeat inference',
                             'conservatism can lose small-motion resolution',
                             'heterogeneous motion and window weights can change the target relative to arithmetic',
                             'no automatic variant selector or physical validation follows'],
                all_frequencies_usable=bool(reports and all(f['usable'] for f in reports)),
                physical_assumptions_verified=False, per_frequency=reports)
