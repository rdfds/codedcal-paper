"""Numerical estimand checks for explicitly supplied linear observation models.

This is standard null-space algebra, not experimental validation. Parameters
and observations must use documented units; near-singular conclusions depend
on the reported relative SVD tolerance. A replicate flag is not a model.
"""
import numpy as np


def assess(observation, estimand, rtol=1e-10):
    """Check ker(A) subset ker(C) for y=A theta and target C theta.

For a failed target return a unit parameter perturbation invisible to A but
visible to C. For a passed target return weights L with L A approximately C.
These are inspectable algebraic witnesses under the specified model only.
"""
    A, C = np.asarray(observation, float), np.asarray(estimand, float)
    if (A.ndim != 2 or C.ndim != 2 or not A.shape[1] or
            C.shape[1] != A.shape[1] or not C.shape[0] or
            not np.isfinite(A).all() or not np.isfinite(C).all()):
        raise ValueError('finite observation and estimand matrices must share columns')
    if not np.isfinite(rtol) or not 0 < rtol < 1:
        raise ValueError('rtol must be finite and between zero and one')
    _, s, vt = np.linalg.svd(A, full_matrices=True)
    cutoff = rtol * s[0] if len(s) else 0.
    rank = int(np.sum(s > cutoff))
    null = vt[rank:].T
    leakage = C @ null
    norms = np.linalg.norm(C, axis=1)
    relative = np.linalg.norm(leakage, axis=1) / np.maximum(norms, 1e-300)
    identified = relative <= rtol
    records = []
    for i, ok in enumerate(identified):
        if ok:
            weights = C[i] @ np.linalg.pinv(A, rcond=rtol)
            records.append(dict(identified=True, observation_weights=weights.tolist(),
                                reconstruction_residual=float(np.linalg.norm(weights @ A-C[i]))))
        else:
            witness = null @ leakage[i]
            witness /= np.linalg.norm(witness)
            records.append(dict(identified=False, parameter_witness=witness.tolist(),
                observation_residual=float(np.linalg.norm(A @ witness)),
                target_change=float(C[i] @ witness)))
    return dict(status='NUMERICAL MODEL CHECK; physical assumptions unverified',
        rank=rank, n_parameters=A.shape[1], nullity=A.shape[1]-rank,
        relative_svd_tolerance=rtol, singular_values=s.tolist(),
        all_targets_identified=bool(np.all(identified)),
        target_null_projection_relative=relative.tolist(), targets=records)


def relative_contrasts(n):
    """Independent contrasts of the first n-1 parameters to the last."""
    if not isinstance(n, (int, np.integer)) or n < 2:
        raise ValueError('at least two parameters required')
    return np.column_stack([np.eye(n-1), -np.ones(n-1)])
