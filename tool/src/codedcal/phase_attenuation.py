"""Phase-preserving matched-count attenuation, with explicit confidence support.

An extension of the existing permanent-code observable, not a new principle of
complex coherence. Under a common displacement characteristic m(q) and a
nonzero source-uniform kernel change psi_sB=c(q)*psi_sA, complex window ratios
obey r_A=m*conj(c), r_B=conj(m)/conj(c). Thus r_A*r_B=|m|**2, even with
unequal window weights. Real projection before multiplication loses this
identity. The target is |m|, not mean displacement, signed attenuation, total
molecular error or a decomposition of common and source-specific motion.

Four component Fieller sets form a Bonferroni rectangle. Marginal coverage is
approximate for independent group means, as in the existing estimator. Full
component sets and unresolved support are retained; no data-selected switching.
"""
import numpy as np

from .smlm import H, PERMS, TRIPLES_A, TRIPLES_B
from .uncertainty import _interval, fieller_set


def phase_components(z, n):
    """One frequency: (groups, codes, panels) -> (groups, windows, Re C/Im C/R).

    Preserves the historical phase/count placements, including both count
    panels and both cross orientations. Only the within-window reference is
    projected to real. Existing real-valued routines are not modified.
    """
    z, n = np.asarray(z), np.asarray(n)
    if z.ndim != 3 or z.shape[1:] != (4, 4) or n.shape != z.shape:
        raise ValueError('z and n must have shape (groups, 4, 4) at one frequency')
    if np.iscomplexobj(n) or not np.isfinite(z).all() or not np.isfinite(n).all() or np.any(n < 0):
        raise ValueError('z must be finite; counts must be finite real and nonnegative')
    a = np.einsum('gct,cb->gtb', z, H)
    b = np.einsum('gct,cb->gtb', n, H)

    def t3(x, y, w):
        return sum(x[:, p[0]]*y[:, p[1]]*w[:, p[2]] for p in PERMS)/6

    result = []
    with np.errstate(over='ignore', invalid='ignore'):
        for triples in (TRIPLES_A, TRIPLES_B):
            reference = 0.; cross = 0.
            for i, j, k in triples:
                reference = reference + t3(a[:, i], a[:, j].conj(), b[:, k]).real/2
                cross = cross + (t3(a[:, i], b[:, j], a[:, k].conj())
                                 + t3(b[:, i], a[:, j], a[:, k].conj()))/4
            result.append(np.stack([cross.real, cross.imag, reference], axis=-1))
    result = np.stack(result, axis=1)
    if not np.isfinite(result).all():
        raise ValueError('nonfinite cubic calculation; rescale all aggregate weights consistently')
    return result


def _merge(parts):
    """Union of nonnegative intervals; preserve open endpoints and gaps."""
    merged = []
    for part in sorted(parts, key=lambda p: p['lower']):
        if not merged:
            merged.append(part.copy()); continue
        last = merged[-1]
        if last['upper'] is not None and (part['lower'] > last['upper'] or
                (part['lower'] == last['upper'] and not
                 (part['lower_closed'] or last['upper_closed']))):
            merged.append(part.copy()); continue
        if last['upper'] is None:
            continue
        if part['upper'] is None or part['upper'] > last['upper']:
            last['upper'], last['upper_closed'] = part['upper'], part['upper_closed']
        elif part['upper'] == last['upper']:
            last['upper_closed'] |= part['upper_closed']
    return merged


def _absolute_bounds(part):
    lo = -np.inf if part['lower'] is None else part['lower']
    hi = np.inf if part['upper'] is None else part['upper']
    if lo < 0 < hi:
        minimum, closed_min = 0., True
    elif lo >= 0:
        minimum, closed_min = lo, part['lower_closed']
    else:
        minimum, closed_min = -hi, part['upper_closed']
    if not np.isfinite([lo, hi]).all():
        return minimum, None, closed_min, False
    if abs(lo) > abs(hi):
        maximum, closed_max = abs(lo), part['lower_closed']
    elif abs(hi) > abs(lo):
        maximum, closed_max = abs(hi), part['upper_closed']
    else:
        maximum, closed_max = abs(lo), part['lower_closed'] or part['upper_closed']
    return minimum, maximum, closed_min, closed_max


def _finite(value):
    if not np.isfinite(value):
        raise ArithmeticError('nonfinite confidence calculation')
    return float(value)


def _modulus_image(real, imag):
    if real is None or imag is None:
        return None
    parts = []
    for x in real:
        for y in imag:
            lx, ux, clx, cux = _absolute_bounds(x)
            ly, uy, cly, cuy = _absolute_bounds(y)
            lower = _finite(np.hypot(lx, ly))
            upper = None if ux is None or uy is None else _finite(np.hypot(ux, uy))
            parts.append(_interval(lower, upper, clx and cly,
                                   upper is not None and cux and cuy))
    return _merge(parts)


def _root_product(left, right):
    if left is None or right is None:
        return None
    parts = []
    for a in left:
        for b in right:
            lo = _finite(np.sqrt(a['lower'])*np.sqrt(b['lower']))
            closed_lo = ((a['lower'] == 0 and a['lower_closed']) or
                         (b['lower'] == 0 and b['lower_closed'])) if lo == 0 else (
                             a['lower_closed'] and b['lower_closed'])
            if a['upper'] == 0 or b['upper'] == 0:
                hi, closed_hi = 0., True
            elif a['upper'] is None or b['upper'] is None:
                hi, closed_hi = None, False
            else:
                hi = _finite(np.sqrt(a['upper'])*np.sqrt(b['upper']))
                closed_hi = a['upper_closed'] and b['upper_closed']
            parts.append(_interval(lo, hi, closed_lo, closed_hi))
    return _merge(parts)


def phase_ratio_set(statistics, confidence_level=.95):
    """Set for sqrt(|(E C_A/E R_A)*(E C_B/E R_B)|) at one q.

    Real input shape: (independent groups, two windows, Re C/Im C/R).
    Negative and zero real/imaginary components are legitimate, not failed
    attenuation estimates. Positive references and strictly positive finite
    magnitude support are required for `usable`; full sets remain visible.
    """
    if not np.isfinite(confidence_level) or not 0 < confidence_level < 1:
        raise ValueError('confidence_level must be finite and between 0 and 1')
    component_level = 1-(1-float(confidence_level))/4
    if component_level >= 1:
        raise ValueError('confidence_level is too close to 1 for four-component allocation')
    if np.iscomplexobj(statistics):
        raise ValueError('statistics must be real columns Re C, Im C, R')
    data = np.asarray(statistics, float)
    if data.ndim != 3 or data.shape[1:] != (2, 3):
        raise ValueError('statistics must have shape (groups, 2, 3)')
    windows = []; numerical_failure = False
    for w in range(2):
        fs = [fieller_set(data[:, w, i], data[:, w, 2], component_level) for i in range(2)]
        for f in fs:
            if any(p[k] is not None and not np.isfinite(p[k])
                   for p in f['intervals'] or [] for k in ('lower', 'upper')):
                f.update(kind='invalid', intervals=None, usable=False,
                         reason='nonfinite_confidence_calculation')
            if f['ratio'] is not None and not np.isfinite(f['ratio']):
                f['ratio'] = None
            # These are Cartesian coordinates, not positive attenuations.
            f['usable'] = bool(f['kind'] == 'bounded' and f['positive_reference'])
            if f['usable']:
                f['reason'] = None
        try:
            magnitude = _modulus_image(fs[0]['intervals'], fs[1]['intervals'])
        except ArithmeticError:
            magnitude = None; numerical_failure = True
        windows.append(dict(real=fs[0], imag=fs[1], modulus_intervals=magnitude))
    try:
        with np.errstate(over='ignore', invalid='ignore'):
            parts = _root_product(*(w['modulus_intervals'] for w in windows))
    except ArithmeticError:
        parts = None; numerical_failure = True
    components = [w[axis] for w in windows for axis in ('real', 'imag')]
    positive_reference = all(f['positive_reference'] for f in components)
    bounded_components = all(f['kind'] == 'bounded' for f in components)
    finite_support = bool(parts and all(p['upper'] is not None for p in parts))
    usable = bool(positive_reference and bounded_components and finite_support and
                  len(parts) == 1 and parts[0]['lower'] > 0)
    ratio = None
    if positive_reference and all(f['ratio'] is not None for f in components):
        moduli = [np.hypot(w['real']['ratio'], w['imag']['ratio']) for w in windows]
        with np.errstate(over='ignore', invalid='ignore'):
            point = np.sqrt(moduli[0])*np.sqrt(moduli[1])
        if np.isfinite(point):
            ratio = float(point)
    contains_one = None if parts is None else any(p['lower'] <= 1 and
                      (p['upper'] is None or p['upper'] >= 1) for p in parts)
    kind = ('invalid' if parts is None else 'empty' if not parts else
            'disjoint' if len(parts) > 1 else
            'bounded' if parts[0]['upper'] is not None else 'unbounded_nonnegative')
    reason = (None if usable else 'nonfinite_confidence_calculation' if numerical_failure
              else 'invalid_component_confidence_set' if any(f['kind'] == 'invalid' for f in components)
              else 'nonpositive_window_reference' if not positive_reference
              else 'unbounded_component_confidence_set' if not bounded_components
              else 'magnitude_not_separated_from_zero' if parts and parts[0]['lower'] == 0
              else 'empty_or_unbounded_magnitude_set')
    return dict(confidence_level=float(confidence_level), component_confidence_level=component_level,
                n_groups=len(data), ratio=ratio, kind=kind, usable=usable, reason=reason,
                positive_reference=positive_reference, intervals=parts, contains_one=contains_one,
                windows=windows)


def phase_attenuation_uncertainty(statistics, q, confidence_level=.95):
    """Pointwise sets for real statistics (groups, frequencies, windows, 3)."""
    if np.iscomplexobj(q) or np.iscomplexobj(statistics):
        raise ValueError('q and statistics must be real')
    q, data = np.asarray(q, float), np.asarray(statistics, float)
    if q.ndim != 1 or not len(q) or not np.isfinite(q).all() or np.any(q <= 0):
        raise ValueError('q must contain finite positive spatial frequencies')
    if data.ndim != 4 or data.shape[1:] != (len(q), 2, 3):
        raise ValueError('statistics must have shape (groups, len(q), 2, 3)')
    rows = []
    for j, frequency in enumerate(q):
        row = phase_ratio_set(data[:, j], confidence_level)
        transformed = []; transform_reason = None

        def gaussian(value):
            log_value = -np.log(value)
            # Sequential division avoids squaring q outside floating support.
            with np.errstate(over='ignore', under='ignore', divide='ignore'):
                result = float((log_value/frequency)/frequency)
            if not np.isfinite(result) or (result == 0 and log_value != 0):
                raise ArithmeticError('nonfinite or underflowed Gaussian transform')
            return result

        for part in row['intervals'] or []:
            lo, hi = part['lower'], part['upper']
            if hi is not None and hi <= 0:
                continue
            try:
                transformed.append(_interval(gaussian(hi) if hi is not None else None,
                    gaussian(lo) if lo > 0 else None,
                    hi is not None and part['upper_closed'], lo > 0 and part['lower_closed']))
            except ArithmeticError:
                transformed = None; transform_reason = 'nonrepresentable_gaussian_transform'; break
        row.update(q=float(frequency), gaussian_equivalent_T_intervals=(
            None if row['kind'] == 'invalid' else transformed),
            gaussian_transform_reason=transform_reason,
            unit_attenuation_excluded=bool(row['usable'] and not row['contains_one']),
            positive_motion_resolved=bool(row['usable'] and row['intervals'][0]['upper'] < 1),
            attenuation_model_compatible=(None if row['kind'] == 'invalid' else any(
                p['lower'] <= 1 and (p['upper'] is None or p['upper'] > 0)
                for p in row['intervals'])))
        rows.append(row)
    return dict(method='phase-preserving matched-count magnitude; four-component Fieller Bonferroni set',
                target='magnitude of the common displacement characteristic under the stated factorization',
                confidence_level=float(confidence_level),
                component_confidence_level=1-(1-float(confidence_level))/4,
                interval_endpoint_convention='null lower/upper means unbounded below/above',
                assumptions=['independent groups; paired windows/components may be dependent',
                    'permanent source-code exogeneity and conditional panel independence',
                    'stable kernels within each window, with nonzero source-uniform complex change between windows',
                    'a common displacement characteristic conditional on the count/kernel weights'],
                limitations=['approximate component coverage; conservative pointwise joint rectangle',
                    'magnitude of a sample ratio is biased; no efficiency or finite-sample unbiasedness claim',
                    'source-dependent changes, code dependence and shared errors remain unresolved',
                    'mean displacement and kernel phase are confounded; neither is recovered',
                    'heterogeneous motion/weights can change the target; internal groups do not resample an acquisition-wide shift',
                    'Gaussian variance interpretation requires a Gaussian displacement law; total molecular accuracy is not identified',
                    'not independently physically validated; no automatic variant selection'],
                all_frequencies_usable=bool(rows and all(r['usable'] for r in rows)),
                physical_assumptions_verified=False, per_frequency=rows)


def phase_budget(z, n, q, confidence_level=.95):
    """Compute the opt-in phase-preserving result from retained complex aggregates."""
    if np.iscomplexobj(q):
        raise ValueError('q must be real')
    q = np.asarray(q, float); z, n = np.asarray(z), np.asarray(n)
    if q.ndim != 1 or not len(q) or not np.isfinite(q).all() or np.any(q <= 0):
        raise ValueError('q must contain finite positive spatial frequencies')
    if z.ndim != 4 or z.shape[1:] != (4, 4, len(q)) or n.shape != z.shape[:3]:
        raise ValueError('z must have shape (groups, 4, 4, len(q)); counts must match')
    if len(z) < 2:
        raise ValueError('at least two independent groups are required for uncertainty')
    comp = np.stack([phase_components(z[..., j], n) for j in range(len(q))], axis=1)
    return dict(n_groups=len(z), q=q.tolist(),
                phase_uncertainty=phase_attenuation_uncertainty(comp, q, confidence_level))
