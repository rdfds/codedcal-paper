"""Bounded image-only observations for a future workflow recommendation.

The summaries do not determine target attainment by themselves. A decision
model requires separate output-recovery training and external evaluation.
Inputs are fixed, nonnegative photometrically normalized post-acquisition crops;
their normalization cost belongs in the pilot budget. No fitted candidate map,
reference calls, injected displacement or outcome is accepted by this API.
"""
import numpy as np
from scipy.ndimage import maximum_filter


def observed_pilot_summary(patches, codebook, minimum_norm=2., max_distance=.6,
                           max_cells=2048, maximum_side=32):
    """Inspect every cell in preselected small crops, with a hard data bound.

    Ambiguity is the gap between the two largest normalized word projections.
    It is not a source-identifiability certificate. Peak counts are an image
    crowding proxy, not a molecular count. Zero signal is unsupported.
    """
    patches = np.asarray(patches, float); book = np.asarray(codebook, float)
    if (patches.ndim != 4 or book.ndim != 2 or len(book) < 2
            or patches.shape[1] != book.shape[1]
            or not np.isfinite(patches).all() or not np.isfinite(book).all()
            or np.any(patches < 0) or np.any(book < 0)
            or np.any(np.linalg.norm(book, axis=1) == 0)
            or not np.isfinite(minimum_norm) or minimum_norm <= 0
            or not np.isfinite(max_distance) or not 0 <= max_distance <= 2
            or not isinstance(max_cells, int) or max_cells < 1
            or not isinstance(maximum_side, int) or maximum_side < 1):
        raise ValueError('finite nonnegative bounded crops and matching codebook required')
    cells = patches.shape[0]*patches.shape[2]*patches.shape[3]
    if (cells < 1 or cells > max_cells or patches.shape[2] > maximum_side
            or patches.shape[3] > maximum_side):
        raise ValueError('pilot exceeds its prespecified image-cell/side budget')
    traces = np.moveaxis(patches, 1, -1).reshape(cells, len(book.T))
    norms = np.linalg.norm(traces, axis=1)
    unit = traces/np.maximum(norms[:, None], np.finfo(float).tiny)
    words = book/np.linalg.norm(book, axis=1)[:, None]
    similarities = unit @ words.T
    top = np.partition(similarities, -2, axis=1)[:, -2:]
    best = top.max(axis=1); gap = best-top.min(axis=1)
    distance = np.sqrt(np.maximum(2-2*best, 0))
    bright = norms >= minimum_norm
    peak_counts = []
    for patch in patches:
        image = np.linalg.norm(patch, axis=0)
        peaks = (image == maximum_filter(image, size=3, mode='constant')) & (image >= minimum_norm)
        # Flat bright plateaus are not counted as many independent peaks.
        from scipy.ndimage import label
        _, count = label(peaks)
        peak_counts.append(int(count))
    summaries = dict(signal_cell_fraction=float(bright.mean()),
                     word_supported_cell_fraction=float((bright & (distance <= max_distance)).mean()),
                     median_bright_word_distance=float(np.median(distance[bright])) if bright.any() else None,
                     median_bright_word_gap=float(np.median(gap[bright])) if bright.any() else None,
                     norm_p95=float(np.quantile(norms, .95)),
                     readout_p95=np.quantile(traces, .95, axis=0).tolist(),
                     image_peaks=peak_counts)
    return dict(supported=bool(bright.any()), summaries=summaries,
                budget=dict(crops=len(patches), cells=cells, measured_values=int(patches.size),
                            candidate_maps_fitted=0, full_field_decoders_run=0),
                input_contract='Fixed post-acquisition photometric crops and codebook only',
                interpretation='Observable pilot features; not a validated workflow recommendation or molecular truth')
