"""Robust readout-response fitting from a separate source calibration region.

Source identities supplied here are calibration proposals, not molecular truth.
Readouts without adequate positive on-source support remain unqualified; their
neutral numerical prior is never a positive source-validation certificate.
"""
import numpy as np


def fit_readout_response(traces, assigned_codes, minimum_sources=3, iterations=100):
    traces = np.asarray(traces, float); codes = np.asarray(assigned_codes, float)
    if (traces.ndim != 2 or traces.shape != codes.shape
            or not np.isfinite(traces).all() or not np.isfinite(codes).all()
            or not np.all((codes == 0) | (codes == 1))
            or not isinstance(minimum_sources, int) or minimum_sources < 2
            or not isinstance(iterations, int) or iterations < 1):
        raise ValueError('finite matched binary calibration traces required')
    positive = (codes > 0) & (traces > 0)
    support = positive.sum(axis=0); qualified = support >= minimum_sources
    # A single amplitude scale needs connected readout/source support. Merely
    # seeing each channel three times does not identify relative gains between
    # disconnected groups of sources.
    active = np.flatnonzero(qualified); groups = []
    adjacency = (positive[:, active].T.astype(int) @ positive[:, active].astype(int)) > 0
    pending = set(range(len(active)))
    while pending:
        todo = [min(pending)]; group = set()
        while todo:
            node = todo.pop()
            if node in group:
                continue
            group.add(node); pending.discard(node)
            todo.extend(set(np.flatnonzero(adjacency[node]))-group)
        groups.append(active[sorted(group)].tolist())
    if len(groups) > 1:
        qualified[:] = False
    gains = np.ones(traces.shape[1]); amplitude = np.zeros(len(traces))
    for iteration in range(1, iterations+1):
        for i in range(len(traces)):
            use = positive[i] & qualified
            amplitude[i] = np.median(traces[i, use]/gains[use]) if use.any() else 0.
        new = np.ones_like(gains)
        for b in np.flatnonzero(qualified):
            use = positive[:, b] & (amplitude > 0)
            if use.sum() >= minimum_sources:
                new[b] = np.median(traces[use, b]/amplitude[use])
            else:
                qualified[b] = False
        if qualified.any():
            new[qualified] /= np.median(new[qualified])
        change = float(np.max(np.abs(new-gains))) if len(gains) else 0.
        gains = new
        if change <= 1e-8:
            break
    # Dispersion is diagnostic, not a calibrated interval or dropout model.
    dispersion = [float(np.median(np.abs(traces[positive[:, b], b]/
                   np.maximum(amplitude[positive[:, b]], 1e-12)-gains[b])))
                  if qualified[b] else None for b in range(len(gains))]
    return dict(gains=gains.tolist(), qualified=qualified.tolist(),
                positive_source_support=support.tolist(), median_absolute_ratio_deviation=dispersion,
                minimum_sources=minimum_sources, iterations=iteration,
                connected_readout_groups=groups, connected_support=len(groups) == 1,
                interpretation='Separate-region robust response diagnostic, not physical gain truth or probability calibration')
