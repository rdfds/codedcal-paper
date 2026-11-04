"""Source checks whose spatial/code fit excludes the tested readouts.

These checks predict measured prepared images. They do not give independent
RNA truth or a calibrated discovery rate. Photometric scales/PSF must already
be fixed outside the checked spatial region.
"""
import numpy as np
from .processed_sparse import solve
from .mixture_quality import multiword_components, tiled_forward


def fit_fold_sources(stack, codebook, operator, heldout, penalty=.5,
                     iterations=150, minimum_norm=2., max_relative_error=.6,
                     border=8):
    """Fit on other readouts, then check a source's withheld contribution.

    The predicted stack, positions and labels never use heldout measurements.
    Validation residuals affect acceptance only. Zero-support or duplicate
    normalized fitting words cannot receive identified-label support.
    """
    stack = np.asarray(stack, float); book = np.asarray(codebook, float)
    held = np.asarray(heldout)
    if (stack.ndim != 3 or book.ndim != 2 or book.shape[1] != len(stack)
            or held.ndim != 1 or not np.issubdtype(held.dtype, np.integer)
            or not len(held) or len(np.unique(held)) != len(held)
            or np.any(held < 0) or np.any(held >= len(stack))
            or len(held) == len(stack) or not np.isfinite(stack).all()
            or not np.isfinite(book).all() or np.any(book < 0)
            or not np.isfinite(minimum_norm) or minimum_norm < 0
            or not np.isfinite(max_relative_error) or max_relative_error < 0):
        raise ValueError('valid finite readout fold and matching codebook required')
    fitting = np.setdiff1d(np.arange(len(stack)), held)
    norms = np.linalg.norm(book[:, fitting], axis=1)
    supported = norms > 0
    unit = book[supported][:, fitting] / norms[supported, None]
    _, groups, counts = np.unique(np.round(unit, 12), axis=0,
                                  return_inverse=True, return_counts=True)
    identified = np.zeros(len(book), bool)
    identified[supported] = counts[groups] == 1
    side = int(np.sqrt(np.asarray(operator).shape[0]))
    if (np.asarray(operator).shape != (side*side, side*side)
            or stack.shape[1] % side or stack.shape[2] % side or not supported.any()):
        raise ValueError('supported pixel-grid fold tiles required')
    density = np.zeros((*stack.shape[1:], len(book)), float)
    fitting_book = book[supported][:, fitting]
    diagnostics = []
    for y in range(0, stack.shape[1], side):
        for x in range(0, stack.shape[2], side):
            measured = np.moveaxis(stack[fitting, y:y+side, x:x+side], 0, -1)
            values, info = solve(operator, fitting_book,
                                 measured.reshape(side*side, len(fitting), order='F'),
                                 penalty=penalty, iterations=iterations)
            tile = values.reshape(side, side, int(supported.sum()), order='F')
            density[y:y+side, x:x+side, supported] = tile
            diagnostics.append(info['proximal_fixed_point_relative'])
    prediction = tiled_forward(density @ book[:, held], operator)
    residual = np.moveaxis(stack[held], 0, -1) - prediction
    components = multiword_components(density, min_area=1, border=border)
    rows = []
    for c in components:
        word = c['word']
        if (not identified[word] or not np.any(book[word, held] > 0)
                or c['peak']*np.linalg.norm(book[word]) < minimum_norm):
            continue
        latent = np.zeros(stack.shape[1]*stack.shape[2])
        latent[c['indices']] = c['coefficients']/c['peak']
        template = tiled_forward(latent.reshape(*stack.shape[1:], 1), operator)[:, :, 0]
        energy = float(np.sum(template**2))
        if energy == 0:
            continue
        support = np.flatnonzero(template)
        projection = template.ravel()[support] @ residual.reshape(-1, len(held))[support]
        expected = c['peak']*book[word, held]
        observed = projection/energy + expected
        relative_error = float(np.linalg.norm(observed-expected)/np.linalg.norm(expected))
        gain = float(energy*(expected@expected) + 2*(projection@expected))
        accepted = gain > 0 and relative_error <= max_relative_error
        rows.append(dict(word=word, position=c['position'], peak=c['peak'],
                         heldout_on_readouts=int((book[word, held] > 0).sum()),
                         expected_trace=expected.tolist(), measured_trace=observed.tolist(),
                         relative_error=relative_error, deletion_sse_gain=gain,
                         accepted=bool(accepted)))
    return dict(fitting_readouts=fitting.tolist(), heldout_readouts=held.tolist(),
                unsupported_words=int((~supported).sum()),
                colliding_fitting_words=int((supported & ~identified).sum()),
                source_rows=rows, prediction=prediction, density=density,
                maximum_fixed_point_relative=max(diagnostics),
                interpretation='Withheld-readout model prediction; not molecule truth or FDR')


def match_crossfit_support(positions, labels, folds, minimum_folds=2, radius=1.5):
    """Count distinct positive withheld folds near each proposed source.

    A distinct fold can contain multiple on-readouts but counts once. This is
    stability evidence; multiple proposed sources cannot reuse a fold source.
    """
    from scipy.optimize import linear_sum_assignment
    xy = np.asarray(positions, float); labels = np.asarray(labels)
    if (xy.shape != (len(labels), 2) or labels.ndim != 1
            or not np.issubdtype(labels.dtype, np.integer) or not np.isfinite(xy).all()
            or not isinstance(minimum_folds, int) or minimum_folds < 1
            or not np.isfinite(radius) or radius <= 0):
        raise ValueError('valid source matching inputs required')
    hits = np.zeros(len(labels), int)
    for fold in folds:
        rows = [r for r in fold['source_rows'] if r['accepted']]
        for word in np.unique(labels):
            ids = np.flatnonzero(labels == word)
            available = [r['position'] for r in rows if r['word'] == word]
            if not available:
                continue
            targets = np.asarray(available)
            # Components for one word are few; avoid a large cross-word matrix.
            distances = np.linalg.norm(xy[ids, None]-targets[None], axis=-1)
            reward = (min(distances.shape)+1)*radius
            rr, cc = linear_sum_assignment(np.where(distances <= radius, distances-reward, reward))
            valid = distances[rr, cc] <= radius
            hits[ids[rr[valid]]] += 1
    return hits >= minimum_folds, hits
