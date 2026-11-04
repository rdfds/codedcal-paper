"""Predictive checks across withheld readouts, not physical molecular truth.

The density/word fit never sees the held-out pixel values. These checks expose
whether a specified code model predicts other observations; correlated noise,
interpolation and model misspecification can still affect their interpretation.
"""
from __future__ import annotations

import numpy as np


def conditional_reconstruction(stack, codebook, heldout):
    """Fit one nonnegative code amplitude per pixel using other readouts only.

    Returns predictions and support diagnostics. Codebook rows are words;
    stack shape is (readouts, height, width). Unsupported and colliding training
    words are explicitly diagnosed rather than treated as identified sources.
    This is an evaluation model, not a replacement for a production decoder.
    """
    stack = np.asarray(stack, float)
    words = np.asarray(codebook, float)
    heldout = np.asarray(heldout, int)
    if stack.ndim != 3 or words.ndim != 2 or words.shape[1] != len(stack):
        raise ValueError("inconsistent stack/codebook dimensions")
    if not np.isfinite(stack).all() or not np.isfinite(words).all():
        raise ValueError("nonfinite prediction input")
    if np.any(stack < 0) or np.any(words < 0):
        raise ValueError("nonnegative observations/codebook required")
    if len(heldout) == 0 or len(np.unique(heldout)) != len(heldout):
        raise ValueError("nonempty unique heldout readouts required")
    if np.any(heldout < 0) or np.any(heldout >= len(stack)):
        raise ValueError("heldout index out of range")
    training = np.setdiff1d(np.arange(len(stack)), heldout)
    if not len(training):
        raise ValueError("no fitting readouts remain")
    matrix = words[:, training]
    norm2 = (matrix ** 2).sum(1)
    available = norm2 > 0
    if not available.any():
        raise ValueError("no supported word remains")
    # Identical normalized traces can represent different amplitudes as well
    # as different words. Diagnose normalized, not just binary, collisions.
    normalized = np.zeros_like(matrix)
    normalized[available] = matrix[available] / np.sqrt(norm2[available, None])
    _, group, counts = np.unique(normalized, axis=0, return_inverse=True,
                                  return_counts=True)
    identifiable = available & (counts[group] == 1)
    x = stack[training].reshape(len(training), -1).T
    similarities = x @ normalized.T
    similarities[:, ~available] = -np.inf
    choice = np.argmax(similarities, axis=1)
    amplitude = np.maximum((x * matrix[choice]).sum(1) / norm2[choice], 0)
    prediction = (words[choice][:, heldout] * amplitude[:, None]).T
    prediction = prediction.reshape((len(heldout),) + stack.shape[1:])
    supported = identifiable[choice].reshape(stack.shape[1:])
    return prediction, dict(training_readouts=training.tolist(),
        heldout_readouts=heldout.tolist(),
        unsupported_words=int((~available).sum()),
        colliding_supported_words=int((available & ~identifiable).sum()),
        identified_pixel_mask=supported,
        fitted_word_indices=choice.reshape(stack.shape[1:]))


def predictive_error(stack, codebook, heldout, mask):
    """Sum held-readout squared error on an externally fixed spatial mask.

    Do not derive the mask from whichever calibration/decoder wins. Return
    raw sums to permit a common denominator across competing geometries.
    """
    predicted, diag = conditional_reconstruction(stack, codebook, heldout)
    mask = np.asarray(mask, bool)
    if mask.shape != np.asarray(stack).shape[1:]:
        raise ValueError("mask shape mismatch")
    usable = mask & diag.pop("identified_pixel_mask")
    diag.pop("fitted_word_indices")
    observed = np.asarray(stack, float)[np.asarray(heldout, int)]
    error = float(((observed - predicted)[:, usable] ** 2).sum())
    energy = float((observed[:, usable] ** 2).sum())
    return dict(squared_error=error, observed_energy=energy,
                selected_pixels=int(mask.sum()), supported_pixels=int(usable.sum()),
                **diag)


def transported_prediction_error(prediction, raw_stack, heldout,
                                 raw_to_aligned, origin, mask, support=None):
    """Score predicted intensities against fixed, original held-readout cells.

    Predictions live in aligned coordinates; transport them back to each
    held-readout camera frame. The target values and spatial mask never change
    with calibration. This prevents treating a smoothed target as new evidence.
    Input stack can use fixed photometric preprocessing, but cannot be warped
    differently per candidate. This still is a predictive check, not truth.
    """
    from .evaluation import warp_readouts
    raw=np.asarray(raw_stack,float);heldout=np.asarray(heldout,int)
    prediction=np.asarray(prediction,float);mask=np.asarray(mask,bool)
    matrices=np.asarray(raw_to_aligned,float)[heldout]
    if prediction.shape!=(len(heldout),)+raw.shape[1:] or mask.shape!=raw.shape[1:]:
        raise ValueError('prediction/target dimensions mismatch')
    inverse=matrices.copy()
    inverse[:,:2]=np.linalg.inv(matrices[:,:2])
    inverse[:,2]=-np.einsum('ki,kij->kj',matrices[:,2],inverse[:,:2])
    back=warp_readouts(prediction,inverse,origin)
    supported=np.ones_like(prediction) if support is None else np.broadcast_to(support,prediction.shape).astype(float)
    raw_support=warp_readouts(supported,inverse,origin)
    usable=mask & np.all(raw_support>=1-1e-6,axis=0)
    target=raw[heldout]
    if np.any(mask & ~usable):
        return dict(squared_error=0.,observed_energy=0.,selected_pixels=int(mask.sum()),supported_pixels=0,
                    partially_supported_pixels=int(usable.sum()),reason='missing_original_prediction_support',
                    target_domain='fixed original camera cells; incomplete candidate excluded')
    return dict(squared_error=float(((target-back)[:,usable]**2).sum()),
                observed_energy=float((target[:,usable]**2).sum()),
                selected_pixels=int(mask.sum()),supported_pixels=int(usable.sum()),
                target_domain='fixed original camera cells; predictions transported back')


def choose_calibration(errors, relative_margin=0.03, minimum_pixels=20,
                       baseline="none", priority=None):
    """Choose from observed predictive losses; unavailable cases abstain.

    This contains no injected displacement, oracle map or call-reference
    endpoint. The margin must be set using development data before validation.
    """
    if not 0 <= relative_margin < 1:
        raise ValueError("invalid relative margin")
    names = list(errors) if priority is None else list(priority)
    if set(names) != set(errors) or len(names) != len(set(names)):
        raise ValueError("priority must include every candidate exactly once")
    good = {name: row for name, row in errors.items()
            if row["supported_pixels"] >= minimum_pixels
            and np.isfinite(row["squared_error"])
            and row["observed_energy"] > 0}
    if baseline not in good:
        return dict(choice=baseline, abstained=True, reason="unsupported_baseline")
    if any(row["supported_pixels"] != good[baseline]["supported_pixels"]
           for row in good.values()):
        return dict(choice=baseline, abstained=True, reason="unequal_prediction_support")
    # Common denominator: do not reward a warp for merely lowering its own
    # observed energy. Geometry can still change interpolation/noise; controls
    # and external outcomes are needed to assess that separate limitation.
    best = min(good, key=lambda name: (good[name]["squared_error"], names.index(name)))
    base_error = good[baseline]["squared_error"]
    if base_error <= 0:
        return dict(choice=baseline, abstained=False, reason="zero_baseline_error")
    improvement = 1 - good[best]["squared_error"] / base_error
    selected = best if improvement >= relative_margin else baseline
    return dict(choice=selected, abstained=False, reason="predictive_comparison",
                relative_improvement=float(improvement))
