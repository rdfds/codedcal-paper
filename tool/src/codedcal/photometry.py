"""Shared sparse intensity-vector decoding, independent of temporal calibration.

Supply all coding readouts and scales estimated without validation outcomes.
Blank-derived error proxies and operating-point selection belong to the caller.
This module supplies neither molecular truth nor a calibrated false-discovery rate.
"""
import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates


def _book(codebook):
    w = np.asarray(codebook, dtype=float)
    if w.ndim != 2 or len(w) < 2 or not np.isin(w, [0, 1]).all():
        raise ValueError('codebook must contain at least two binary word rows')
    weight = w.sum(axis=1)
    if weight[0] <= 0 or not np.all(weight == weight[0]):
        raise ValueError('codebook words must have equal positive weight')
    if len(np.unique(w, axis=0)) != len(w):
        raise ValueError('codebook words must be unique')
    return w


def fit_scales(traces, seed_words, codebook, training_mask, min_support=10):
    """Median positive on-bit training signal; stop on unsupported readouts.

    training_mask must exclude validation seeds and include only trusted exact
    words. No inferred word from the scored validation set is used here.
    """
    w = _book(codebook)
    x = np.asarray(traces, dtype=float)
    ids = np.asarray(seed_words)
    train = np.asarray(training_mask)
    if x.ndim != 2 or x.shape[1] != w.shape[1] or not np.isfinite(x).all() or np.any(x < 0):
        raise ValueError('traces must be finite nonnegative seed-by-readout values')
    if ids.shape != (len(x),) or not np.issubdtype(ids.dtype, np.integer) or np.any((ids < 0) | (ids >= len(w))):
        raise ValueError('seed_words must be valid integer codebook row indices')
    if train.shape != (len(x),) or train.dtype != bool:
        raise ValueError('training_mask must be a boolean vector matching seeds')
    if not isinstance(min_support, (int, np.integer)) or min_support < 1:
        raise ValueError('min_support must be a positive integer')
    scales, support = [], []
    for k in range(x.shape[1]):
        use = train & (w[ids, k] > 0) & (x[:, k] > 0)
        support.append(int(use.sum()))
        if support[-1] < min_support:
            raise ValueError('unsupported scale for readout %d: %d positive training seeds' % (k, support[-1]))
        scales.append(float(np.median(x[use, k])))
    return np.array(scales), np.array(support)


def score_traces(traces, codebook, scales):
    """Equal-prior nearest normalized word; equidistant ties abstain."""
    w = _book(codebook)
    x = np.asarray(traces, dtype=float)
    scales = np.asarray(scales, dtype=float)
    if x.ndim != 2 or x.shape[1] != w.shape[1] or not np.isfinite(x).all() or np.any(x < 0):
        raise ValueError('traces must be finite nonnegative seed-by-readout values')
    if scales.shape != (w.shape[1],) or not np.isfinite(scales).all() or np.any(scales <= 0):
        raise ValueError('scales must be finite positive values for every readout')
    v = x / scales
    brightness = np.linalg.norm(v, axis=1)
    if not np.isfinite(brightness).all():
        raise ValueError('scaled traces overflow')
    unit = v / np.maximum(brightness[:, None], 1e-30)
    normalized_words = w / np.sqrt(w.sum(axis=1))[:, None]
    sim = unit @ normalized_words.T
    word = np.argmax(sim, axis=1)
    best = sim[np.arange(len(word)), word]
    runner = np.partition(sim, -2, axis=1)[:, -2]
    return dict(word=word, distance=np.sqrt(np.maximum(0, 2 - 2 * best)),
                brightness=brightness, tie=np.abs(best - runner) <= 1e-10)


def select_calls(centers, scores, source_ids, n_words, distance_cut,
                 brightness_cut, eligible=None, suppression_radius=1.):
    """Select before counting, suppressing nearby same-word duplicates.

    Priority is distance, decreasing brightness, then unique source ID.
    Gene/blank labels never enter selection. Scores must come from score_traces.
    """
    xy = np.asarray(centers, dtype=float)
    ids = np.asarray(source_ids)
    word = np.asarray(scores['word'])
    d, b, tie = (np.asarray(scores[k]) for k in ('distance', 'brightness', 'tie'))
    n = len(word)
    if xy.shape != (n, 2) or not np.isfinite(xy).all():
        raise ValueError('centers must be finite seed-by-two coordinates')
    if ids.shape != (n,) or not np.issubdtype(ids.dtype, np.integer) or len(np.unique(ids)) != n:
        raise ValueError('source_ids must be unique integer identifiers')
    if not isinstance(n_words, (int, np.integer)) or n_words < 2 or word.shape != (n,) or not np.issubdtype(word.dtype, np.integer) or np.any((word < 0) | (word >= n_words)):
        raise ValueError('word indices and n_words do not match')
    if any(v.shape != (n,) for v in (d, b, tie)) or tie.dtype != bool or not np.isfinite(d).all() or not np.isfinite(b).all() or np.any(d < 0) or np.any(b < 0):
        raise ValueError('scores must be finite matching vectors and boolean ties')
    if not np.isfinite([distance_cut, brightness_cut, suppression_radius]).all() or min(distance_cut, brightness_cut, suppression_radius) < 0:
        raise ValueError('cuts and suppression radius must be finite nonnegative values')
    if eligible is None:
        eligible = np.ones(n, dtype=bool)
    eligible = np.asarray(eligible)
    if eligible.shape != (n,) or eligible.dtype != bool:
        raise ValueError('eligible must be a boolean vector matching seeds')
    candidates = np.flatnonzero(eligible & ~tie & (d <= distance_cut) & (b >= brightness_cut))
    order = sorted(map(int, candidates), key=lambda i: (float(d[i]), -float(b[i]), int(ids[i])))
    selected, byword = [], {}
    for i in order:
        others = byword.setdefault(int(word[i]), [])
        if any(np.linalg.norm(xy[i] - xy[j]) <= suppression_radius for j in others):
            continue
        selected.append(i)
        others.append(i)
    selected = np.array(selected, dtype=int)
    return selected, np.bincount(word[selected], minlength=n_words)


def sample_traces(images, centers, raw_to_aligned, origin=(0., 0.),
                  highpass_sigma=2., lowpass_sigma=.6, edge_margin=1.):
    """Sample each readout through its inverse affine map, in supplied order.

    images may be a generator: only one image is preprocessed at a time.
    raw_to_aligned has shape(readouts,3,2): raw @ coef[:2] + coef[2].
    origin is the raw (x,y) origin of each identically sized image crop.
    Return traces and a mask requiring interior sampling in EVERY readout.
    """
    xy = np.asarray(centers, dtype=float)
    coef = np.asarray(raw_to_aligned, dtype=float)
    origin = np.asarray(origin, dtype=float)
    if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
        raise ValueError('centers must be finite seed-by-two coordinates')
    if coef.ndim != 3 or coef.shape[1:] != (3, 2) or not len(coef) or not np.isfinite(coef).all():
        raise ValueError('raw_to_aligned must contain finite readout-by-3-by-2 maps')
    if origin.shape != (2,) or not np.isfinite(origin).all() or not np.isfinite([highpass_sigma, lowpass_sigma, edge_margin]).all() or min(highpass_sigma, lowpass_sigma, edge_margin) < 0:
        raise ValueError('origin, sigmas and margin must be finite; sigmas and margin nonnegative')
    traces = np.zeros((len(xy), len(coef)), dtype=np.float32)
    good = np.ones(len(xy), dtype=bool)
    count, shape = 0, None
    for k, image in enumerate(images):
        if k >= len(coef):
            raise ValueError('more images than affine maps')
        image = np.asarray(image, dtype=np.float32)
        if image.ndim != 2 or not np.isfinite(image).all():
            raise ValueError('images must be finite two-dimensional arrays')
        if shape is not None and image.shape != shape:
            raise ValueError('all readout crops must have identical shape')
        shape = image.shape
        try:
            points = (xy - coef[k, 2]) @ np.linalg.inv(coef[k, :2]) - origin
        except np.linalg.LinAlgError as exc:
            raise ValueError('singular readout affine map') from exc
        high = np.maximum(image - gaussian_filter(image, highpass_sigma), 0)
        smooth = gaussian_filter(high, lowpass_sigma)
        traces[:, k] = map_coordinates(smooth, [points[:, 1], points[:, 0]], order=1, mode='constant', cval=0)
        good &= (points[:, 0] >= edge_margin) & (points[:, 0] <= shape[1]-1-edge_margin) & (points[:, 1] >= edge_margin) & (points[:, 1] <= shape[0]-1-edge_margin)
        count += 1
    if count != len(coef):
        raise ValueError('fewer images than affine maps')
    return traces, good
