"""Require complete source prediction evidence, including off-bit failures."""
import numpy as np
from scipy.optimize import linear_sum_assignment


def complete_singleton_support(positions, labels, folds, readouts, radius=1.5,
                               max_relative_error=.6):
    """Aggregate independently predicted singleton traces without cherry-picking.

    Every readout needs an identified, unambiguous, converged fitting-only
    source at the same site. Both on and off residuals enter the single vector
    error. Missing or unsupported folds cause abstention. This is conditional
    model consistency, not an identity certificate or a discovery rate.
    """
    xy = np.asarray(positions, float); labels = np.asarray(labels)
    if (xy.shape != (len(labels), 2) or not np.isfinite(xy).all()
            or labels.ndim != 1 or not np.issubdtype(labels.dtype, np.integer)
            or not isinstance(readouts, int) or readouts < 2
            or radius <= 0 or not np.isfinite(radius)
            or max_relative_error < 0 or not np.isfinite(max_relative_error)):
        raise ValueError('valid singleton evidence request required')
    expected = np.zeros((len(labels), readouts)); observed = np.zeros_like(expected)
    available = np.zeros_like(expected, bool); seen = set()
    for fold in folds:
        held = fold['heldout_readouts']
        if len(held) != 1 or held[0] in seen or not 0 <= held[0] < readouts:
            raise ValueError('distinct singleton folds required')
        bit = held[0]; seen.add(bit)
        if not fold['selection']['converged'] or not fold['debias']['converged']:
            continue
        rows = [r for r in fold['source_rows'] if r['fitting_identified']
                and not r['ambiguous'] and r['response_supported'] and r.get('profile_supported',True)]
        for word in np.unique(labels):
            ids = np.flatnonzero(labels == word)
            rr = [r for r in rows if r['word'] == word]
            if not rr:
                continue
            target = np.array([r['position'] for r in rr])
            distances = np.linalg.norm(xy[ids, None]-target[None], axis=-1)
            reward = (min(distances.shape)+1)*radius
            a, b = linear_sum_assignment(np.where(distances <= radius, distances-reward, reward))
            for i, j in zip(a, b):
                if distances[i, j] <= radius:
                    source = rr[j]
                    expected[ids[i], bit] = source['expected_trace'][0]
                    observed[ids[i], bit] = source['measured_trace'][0]
                    available[ids[i], bit] = True
    norm = np.linalg.norm(expected, axis=1)
    error = np.linalg.norm(observed-expected, axis=1)/np.maximum(norm, 1e-12)
    keep = available.all(axis=1) & (norm > 0) & (error <= max_relative_error)
    return keep, dict(fold_support=available.sum(axis=1), relative_error=error,
                     expected_trace=expected, measured_trace=observed,
                     complete=available.all(axis=1))
