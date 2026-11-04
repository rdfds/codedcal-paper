"""Pointwise approximate Fieller sets for paired independent-group means.

The inequality is the one used in Source Data 70, validation.py::fieller.
This implementation additionally represents every degenerate set explicitly.
It does not validate physical code exogeneity, source identities or independence.
"""
import numpy as np
from scipy.stats import t as student_t


def _interval(lower=None, upper=None, lower_closed=False, upper_closed=False):
    # None denotes an unbounded endpoint, never a missing finite bound.
    return dict(lower=lower, upper=upper,
                lower_closed=lower_closed, upper_closed=upper_closed)


def _quadratic_set(a, b, c):
    """Solve a*r**2 + b*r + c <= 0, without imposing r >= 0."""
    if a == 0:
        if b == 0:
            return ('all_real', [_interval()]) if c <= 0 else ('empty', [])
        root = float(-c / b)
        return 'half_line', ([_interval(upper=root, upper_closed=True)] if b > 0
                             else [_interval(lower=root, lower_closed=True)])
    disc = b*b - 4*a*c
    # Resolve cancellation at a double root, rather than turning a point set
    # into an empty set through a one-ulp negative discriminant.
    if abs(disc) <= 8*np.finfo(float).eps*(b*b + abs(4*a*c)):
        disc = 0.
    if disc < 0:
        return ('all_real', [_interval()]) if a < 0 else ('empty', [])
    if disc == 0:
        if a < 0:
            return 'all_real', [_interval()]
        root = float(-b / (2*a))
        return 'bounded', [_interval(root, root, True, True)]
    # Stable roots when one of -b +/- sqrt(disc) suffers cancellation.
    v = -.5 * (b + np.copysign(np.sqrt(disc), b))
    lower, upper = sorted((float(v/a), float(c/v)))
    if a > 0:
        return 'bounded', [_interval(lower, upper, True, True)]
    return 'disjoint', [_interval(upper=lower, upper_closed=True),
                        _interval(lower=upper, lower_closed=True)]


def fieller_set(numerator, reference, confidence_level=.95):
    """Approximate confidence set for E[numerator] / E[reference].

    Rows must be independent groups, not individual correlated localizations.
    A Student critical value and the estimated covariance of group means are
    used. Coverage can fail with concentrated weights or model violations.
    Sets are not clipped to positive attenuations or to ratios <= 1.
    """
    if not np.isfinite(confidence_level) or not 0 < confidence_level < 1:
        raise ValueError('confidence_level must be finite and between 0 and 1')
    x, y = np.asarray(numerator, float), np.asarray(reference, float)
    if x.ndim != 1 or y.shape != x.shape:
        raise ValueError('numerator and reference must be matching 1D group arrays')
    n = len(x)
    result = dict(confidence_level=float(confidence_level), n_groups=n,
                  kind='invalid', intervals=None, ratio=None,
                  numerator=None, reference=None, positive_reference=False,
                  usable=False, contains_one=None, reason=None)
    if n < 2 or not (np.isfinite(x).all() and np.isfinite(y).all()):
        result['reason'] = ('fewer_than_two_groups' if n < 2
                            else 'nonfinite_group_statistics')
        return result
    # A common scale leaves the inequality unchanged and prevents count-cubic
    # statistics from overflowing when forming squared means/covariances.
    scale = max(float(np.max(np.abs(x))), float(np.max(np.abs(y)))) or 1.
    values = np.stack([x/scale, y/scale])
    mx, my = values.mean(axis=1)
    cov = np.cov(values, ddof=1)/n
    critical = student_t.ppf((1+confidence_level)/2, n-1)**2
    a = float(my*my-critical*cov[1, 1])
    b = float(-2*mx*my+2*critical*cov[0, 1])
    c = float(mx*mx-critical*cov[0, 0])
    if not np.isfinite([a, b, c]).all():
        result['reason'] = 'nonfinite_confidence_calculation'
        return result
    kind, intervals = _quadratic_set(a, b, c)
    positive_reference = bool(my > 0)
    usable = bool(kind == 'bounded' and intervals[0]['lower'] > 0
                  and positive_reference)
    contains_one = any((i['lower'] is None or i['lower'] <= 1)
                       and (i['upper'] is None or i['upper'] >= 1)
                       for i in intervals)
    result.update(kind=kind, intervals=intervals,
                  ratio=float(mx/my) if my != 0 else None,
                  numerator=float(mx*scale), reference=float(my*scale),
                  positive_reference=positive_reference, usable=usable,
                  contains_one=contains_one,
                  reason=(None if usable else 'nonpositive_reference' if not positive_reference
                          else 'nonpositive_ratio_support' if kind == 'bounded'
                          else kind+'_confidence_set'))
    # Observable group concentration only: these are not the inaccessible
    # source-kernel effective counts in the validation study, nor calibrated
    # independent sample sizes. Signed contributions may cancel in the mean.
    concentration = {}
    for label, contribution in [('numerator', values[0]), ('reference', values[1])]:
        mass = float(np.abs(contribution).sum())
        concentration[label] = dict(
            largest_absolute_fraction=float(np.abs(contribution).max()/mass) if mass else None,
            effective_absolute_contributions=float(mass**2/np.sum(contribution**2)) if mass else None)
    result['group_concentration'] = concentration
    return result


def attenuation_uncertainty(group_statistics, q, confidence_level=.95):
    """Fieller sets and conditional -log(r)/q**2 sets at each frequency.

    Input shape is (groups, frequencies, 2), ordered numerator/reference.
    The transformed set uses only r > 0; the full, possibly negative ratio
    set remains intact. No Gaussian interpretation is imposed on the inputs.
    """
    a, q = np.asarray(group_statistics, float), np.asarray(q, float)
    if q.ndim != 1 or not len(q) or not np.isfinite(q).all() or np.any(q <= 0):
        raise ValueError('q must contain finite positive spatial frequencies')
    if a.ndim != 3 or a.shape[1:] != (len(q), 2):
        raise ValueError('group_statistics must have shape (groups, len(q), 2)')
    reports = []
    for j, frequency in enumerate(q):
        f = fieller_set(a[:, j, 0], a[:, j, 1], confidence_level)
        transformed = []
        for interval in f['intervals'] or []:
            lo, hi = interval['lower'], interval['upper']
            if hi is not None and hi <= 0:
                continue
            transformed.append(_interval(
                lower=float(-np.log(hi)/frequency**2) if hi is not None else None,
                upper=float(-np.log(lo)/frequency**2) if lo is not None and lo > 0 else None,
                lower_closed=interval['upper_closed'] if hi is not None else False,
                upper_closed=interval['lower_closed'] if lo is not None and lo > 0 else False))
        f.update(q=float(frequency), gaussian_equivalent_T_intervals=(
                     None if f['kind'] == 'invalid' else transformed),
                 unit_attenuation_excluded=bool(f['usable'] and not f['contains_one']),
                 positive_motion_resolved=bool(f['usable'] and f['intervals'][0]['upper'] < 1),
                 attenuation_model_compatible=(None if f['kind'] == 'invalid' else any(
                     (i['upper'] is None or i['upper'] > 0)
                     and (i['lower'] is None or i['lower'] <= 1) for i in f['intervals'])))
        reports.append(f)
    return dict(method='pointwise approximate Fieller, Student critical value',
                confidence_level=float(confidence_level),
                interval_endpoint_convention='null lower/upper means unbounded below/above',
                target='ratio of expected matched-count cross and reference moments',
                interpretation='displacement attenuation under the physical model; '
                               '-log(ratio)/q^2 is a Gaussian-equivalent variance only conditionally',
                assumptions=['independent groups', 'physical source-code exogeneity',
                             'conditionally independent temporal panels', 'stable source kernels',
                             'common displacement law for a common variance interpretation'],
                limitations=['approximate pointwise coverage, not simultaneous or exact',
                             'concentrated group contributions can defeat the approximation',
                             'a finite set or quiet D does not validate physical assumptions',
                             'source heterogeneity changes the weighted attenuation target'],
                all_frequencies_usable=bool(reports and all(f['usable'] for f in reports)),
                physical_assumptions_verified=False, per_frequency=reports)
