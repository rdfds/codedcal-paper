"""Memory-bounded evaluation of the retained integer group comparator.

Offsets, codewords, floating-point dtypes and search/tie rules are unchanged.
The original adapter remains available for historical reproduction.
"""
from __future__ import annotations
import numpy as np
from .workflow_validation import warp


def projected_losses(normalized, shifted, mask, inner, offsets, denominator, code_block=32, candidate_batch=9):
    """Evaluate the same full-field maxima with a bounded codeword working set.

    Keep each candidate's complete float32 maximum image. The final float64
    square/sum is therefore the original full-field reduction, not a sum of
    independently rounded tile losses. No codewords or candidate pixels are cut.
    """
    if np.__version__ != '1.26.4':
        raise ValueError('Exact blocked projection requires the qualified NumPy 1.26.4 runtime')
    if (not isinstance(code_block, int) or isinstance(code_block, bool)
            or code_block < 32 or code_block % 32):
        raise ValueError('code_block must be a positive multiple of 32 to preserve tested BLAS blocking')
    if not isinstance(candidate_batch, int) or isinstance(candidate_batch, bool) or candidate_batch < 1:
        raise ValueError('candidate_batch must be a positive integer')
    offsets = list(offsets)
    if not offsets or normalized.ndim != 2 or shifted.ndim != 3:
        raise ValueError('Nonempty offsets and frame/codeword/image arrays required')
    if normalized.shape[0] != shifted.shape[0] or mask.shape != (len(shifted),):
        raise ValueError('Readout dimensions disagree')
    # Reuse these readout views across codeword blocks rather than repeatedly
    # materializing masked full fields at the projection boundary.
    positive = shifted[mask].reshape(int(mask.sum()), -1)
    negative = shifted[~mask].reshape(int((~mask).sum()), -1)
    losses = []
    for first in range(0, len(offsets), candidate_batch):
        current = offsets[first:first+candidate_batch]
        maxima = np.full((len(current),) + tuple(inner), -np.inf, np.float32)
        for start in range(0, normalized.shape[1], code_block):
            stop = min(start + code_block, normalized.shape[1])
            plus = (normalized[mask, start:stop].T @ positive).reshape((-1,) + shifted.shape[1:])
            minus = (normalized[~mask, start:stop].T @ negative).reshape((-1,) + shifted.shape[1:])
            buffer = np.empty((stop - start,) + tuple(inner), np.float32)
            for k, (y, x) in enumerate(current):
                np.add(plus[:, 10-y:10-y+inner[0], 10-x:10-x+inner[1]],
                       minus[:, 10+y:10+y+inner[0], 10+x:10+x+inner[1]], out=buffer)
                np.maximum(maxima[k], buffer.max(axis=0), out=maxima[k])
            # Avoid retaining the previous block during the next allocation.
            del plus, minus, buffer
        losses.extend(1 - float(np.square(image.astype(float)).sum()) / denominator for image in maxima)
        del maxima
    return losses


def striped_native_loss(images, codebook, maps, stripe_rows=16):
    """Evaluate native translation/einsum kernels without a full dense dot tensor.

    Integer maps use the same native Hermite translation operator. Each stripe
    gets its full original-image coordinates (including zero padding), and the
    final native TensorFlow reduction uses the complete maximum image.
    """
    if not isinstance(stripe_rows, int) or isinstance(stripe_rows, bool) or stripe_rows < 1:
        raise ValueError('stripe_rows must be a positive integer')
    import tensorflow as tf
    from bardensr.registration.lowrankregistration import _procminicode
    from bardensr.registration import translations_tf
    if not np.array_equal(maps, np.rint(maps)):
        raise ValueError('The retained group comparator requires integer maps')
    original = images.astype(float)
    mini, code, _, size = _procminicode(original, codebook, 20)
    shape = tuple(int(x) for x in size.numpy())
    maximum = np.empty(shape, np.float64)
    for row in range(0, shape[0], stripe_rows):
        height = min(stripe_rows, shape[0] - row)
        t = tf.convert_to_tensor(-20 - maps + np.array([row, 0]), dtype=tf.float64)
        moved = translations_tf.floating_slices(mini, t, tf.constant([height, shape[1]]), 'hermite')
        dots = tf.einsum('fj,f...->...j', code, moved)
        maximum[row:row+height] = tf.reduce_max(dots, axis=-1).numpy()
        del moved, dots
    # This is the upstream _calc_scaled_loss arithmetic and complete shape.
    loss = -tf.convert_to_tensor(maximum) ** 2
    nloss = tf.cast(tf.reduce_prod(tf.shape(loss)), tf.float64)
    ninput = tf.cast(tf.reduce_prod(tf.shape(mini)), tf.float64)
    return float((1 + (tf.reduce_mean(loss) / tf.reduce_mean(mini ** 2)) * (nloss / ninput)).numpy())


class BoundedGroupMixin:
    """Drop-in group_fit override; all other calibration/decoder code is inherited."""

    def group_fit(self, case, common):
        region = tuple(slice(*r) for r in case.metadata.get('calibration_region',
                       [[0, case.images.shape[1]], [0, case.images.shape[2]], [0, case.images.shape[3]]]))
        images = case.images[(slice(None),) + region]
        alpha = case.support[(slice(None),) + region]
        images = np.where(alpha >= 1-1e-6, images, 0).max(axis=1).astype(np.float32)
        colors = np.unique(case.colors)
        if len(colors) not in (2, 4):
            raise ValueError('Group-search comparator qualified for two/four colors only')
        partitions = [(colors[0],)] if len(colors) == 2 else [(colors[0], c) for c in colors[1:]]
        normalized = (case.codebook / np.sqrt(np.square(case.codebook).sum(axis=0))).astype(np.float32)
        denominator = max(float(np.square(images.astype(float)).sum()), 1e-20)
        best = None
        fits = []
        ysize, xsize = images.shape[1:]
        inner = (ysize+40, xsize+40)
        for start_name, initial in [('zero', np.zeros_like(common[:, 1:])), ('round', np.rint(common[:, 1:]))]:
            initial = initial.astype(int)
            if np.max(np.abs(initial)) > 10:
                continue
            shifted = warp(np.pad(images, ((0, 0), (30, 30), (30, 30)))[:, None],
                           np.c_[np.zeros(len(images)), initial])[:, 0]
            localbest = None
            for partition in partitions:
                mask = np.isin(case.colors, partition)
                def evaluate(offsets):
                    losses = projected_losses(normalized, shifted, mask, inner, offsets, denominator)
                    return [(loss, y*y+x*x, tuple(partition), int(y), int(x))
                            for loss, (y, x) in zip(losses, offsets)]
                found = min(evaluate([(y, x) for y in range(-8, 9, 2) for x in range(-8, 9, 2)]))
                y, x = found[-2:]
                found = min([found] + evaluate([(y+dy, x+dx) for dy in [-1, 0, 1] for dx in [-1, 0, 1]]))
                if localbest is None or found < localbest:
                    localbest = found
            loss, _, partition, y, x = localbest
            sign = np.where(np.isin(case.colors, partition), 1, -1)
            maps = initial + sign[:, None] * np.array([y, x])
            actual = striped_native_loss(images, case.codebook, maps)
            if not np.isfinite(actual) or abs(actual-loss) >= 2e-6:
                raise ValueError('Native group-loss equivalence failed')
            fit = (loss, float(np.square(maps).sum()), start_name, maps, actual)
            fits.append(dict(start=start_name, projected_loss=loss, native_loss=actual))
            if best is None or fit[:3] < best[:3]:
                best = fit
        if best is None:
            raise ValueError('No bounded group-search start')
        return np.c_[np.zeros(len(images)), best[3]], dict(
            calibration_fit_source='two_start_balanced_color_group_search', integer_only_policy=True,
            group_loss=best[4], native_loss_equivalence_checked=True, group_families=fits,
            memory_bounded_group_evaluation=True, code_block=32, candidate_batch=9, native_loss_stripe_rows=16)
