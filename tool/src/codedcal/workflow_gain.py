"""Explicit amplitude gauge for relative gains fitted with cosine objectives.

A cosine-only objective does not identify a common positive gain multiplier.
An absolute intensity cutoff after gain correction therefore requires a declared
reference scale. This helper anchors that otherwise free multiplier; it does
not estimate physical gain, regularize relative gains or qualify a decoder.
"""
from __future__ import annotations
import math
import numpy as np


def _positive_vector(values, name):
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 1 or not result.size or not np.isfinite(result).all() or np.any(result <= 0):
        raise ValueError(name+' must be a nonempty positive finite vector')
    return result


def anchor_relative_gains(gains, reference_gains):
    """Set mean(log(anchored/reference)) to zero without altering gain ratios.

    The reference is mandatory: its amplitude units determine the meaning of a
    downstream absolute threshold. Caller-provided units must be selected before
    evaluating outcomes. No clipping or numerical repair of invalid gains occurs.
    Anchoring leaves cosine distances invariant in exact arithmetic, but changes
    the unanchored absolute-intensity convention. Original results must therefore
    remain original results; this is an explicit alternate inference contract.
    """
    gains = _positive_vector(gains, 'gains')
    reference = _positive_vector(reference_gains, 'reference_gains')
    if gains.shape != reference.shape:
        raise ValueError('gains and reference_gains must have identical shapes')
    log_gains = np.log(gains)
    log_reference = np.log(reference)
    offset = math.fsum((log_gains-log_reference).tolist())/gains.size
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        anchored = np.exp(log_gains-offset)
    if not np.isfinite(anchored).all() or np.any(anchored <= 0):
        raise ValueError('Anchored gains cannot be represented without clipping')
    try:
        scale = math.exp(offset)
    except OverflowError:
        scale = None
    if scale is not None and (not math.isfinite(scale) or scale <= 0):
        scale = None
    diagnostics = dict(
        contract='geometric_mean_gain_ratio_to_declared_reference_equals_one',
        log_common_gain_multiplier=offset,
        common_gain_multiplier=scale,
        log_relative_gain_range=float(np.ptp(log_gains)),
        physical_gain_identified=False,
        relative_gain_stability_established=False,
        downstream_absolute_intensity_units='declared_reference_gains',
    )
    return anchored, diagnostics
