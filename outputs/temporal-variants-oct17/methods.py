"""Specified arithmetic Fieller and geometric Bonferroni-rectangle methods."""
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tool/src'))
from codedcal.uncertainty import fieller_set


def interval(lower, upper, lower_closed=False, upper_closed=False):
    return dict(lower=lower, upper=upper, lower_closed=lower_closed, upper_closed=upper_closed)


def positive_parts(intervals):
    if intervals is None:
        return None
    result = []
    for part in intervals:
        lo, hi = part['lower'], part['upper']
        if hi is not None and hi <= 0:
            continue
        result.append(interval(max(0., lo) if lo is not None else 0., hi,
                               part['lower_closed'] if lo is not None and lo > 0 else False,
                               part['upper_closed']))
    return result


def product_root(left, right):
    """Image of positive interval unions under sqrt(x*y), with no clipping."""
    if left is None or right is None:
        return None
    result = []
    for a in left:
        for b in right:
            lo = float(np.sqrt(a['lower'])*np.sqrt(b['lower']))
            hi = (None if a['upper'] is None or b['upper'] is None else
                  float(np.sqrt(a['upper'])*np.sqrt(b['upper'])))
            result.append(interval(lo, hi, a['lower_closed'] and b['lower_closed'],
                                   hi is not None and a['upper_closed'] and b['upper_closed']))
    result.sort(key=lambda r:r['lower'])
    merged = []
    for part in result:
        if not merged:
            merged.append(part.copy()); continue
        last = merged[-1]
        touches = (last['upper'] is None or part['lower'] < last['upper'] or
                   (part['lower'] == last['upper'] and (last['upper_closed'] or part['lower_closed'])))
        if not touches:
            merged.append(part.copy()); continue
        if last['upper'] is None:
            continue
        if part['upper'] is None or part['upper'] > last['upper']:
            last['upper'], last['upper_closed'] = part['upper'], part['upper_closed']
        elif part['upper'] == last['upper']:
            last['upper_closed'] |= part['upper_closed']
    return merged


def contains(parts, value):
    if parts is None or not np.isfinite(value):
        return False
    for part in parts:
        lo, hi = part['lower'], part['upper']
        tolerance = 8*np.finfo(float).eps*max(1., abs(value), abs(lo or 0.), abs(hi or 0.))
        lower_ok = lo is None or (value >= lo-tolerance if part['lower_closed'] else value > lo)
        upper_ok = hi is None or (value <= hi+tolerance if part['upper_closed'] else value < hi)
        if lower_ok and upper_ok:
            return True
    return False


def geometric_target(population):
    ratios = population[:, 0]/population[:, 1]
    return float(np.exp(np.log(ratios).mean()))


def arithmetic_target(population):
    return float(population[:, 0].sum()/population[:, 1].sum())


def estimate(groups, variant):
    """Only group/window numerator/reference observations enter the estimate."""
    if variant == 'arithmetic':
        pooled = groups.mean(axis=1)
        f = fieller_set(pooled[:, 0], pooled[:, 1])
        return dict(ratio=f['ratio'] if f['positive_reference'] else None,
                    usable=f['usable'], kind=f['kind'], intervals=f['intervals'],
                    windows=None, positive_reference=f['positive_reference'])
    if variant != 'geometric':
        raise ValueError(variant)
    windows = [fieller_set(groups[:, w, 0], groups[:, w, 1], confidence_level=.975)
               for w in range(2)]
    parts = product_root(*(positive_parts(f['intervals']) for f in windows))
    positive_reference = all(f['positive_reference'] for f in windows)
    available = positive_reference and all(f['ratio'] is not None and f['ratio'] > 0 for f in windows)
    ratio = float(np.exp(np.mean([np.log(f['ratio']) for f in windows]))) if available else None
    usable = all(f['usable'] for f in windows)
    kind = ('invalid' if parts is None else 'empty' if not parts else 'disjoint' if len(parts)>1
            else 'bounded' if parts[0]['upper'] is not None else 'unbounded_positive')
    return dict(ratio=ratio, usable=usable, kind=kind, intervals=parts,
                positive_reference=positive_reference,
                windows=[dict(kind=f['kind'], usable=f['usable'], ratio=f['ratio'],
                              positive_reference=f['positive_reference'], intervals=f['intervals']) for f in windows])
