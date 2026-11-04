"""Unpenalized nonnegative amplitude refit after sparse source selection.

Selection and refitting use fitting readouts only. Withheld measurements only
score predictions. Sources/labels/spatial templates remain fixed. This removes
penalty shrinkage from the prediction; it does not establish unique identities.
"""
import numpy as np
from scipy.sparse import csr_matrix, vstack
from .crossfit_quality import fit_fold_sources
from .mixture_quality import multiword_components, tiled_forward


def spectral_gram(templates, codes, block=65536):
    """Sparse separable Gram with bounded temporary spectral allocations.

    Never allocate edges-by-all-readouts for the entire spatial Gram. Blocks
    preserve the same per-edge arithmetic as the original implementation.
    """
    gram = (templates @ templates.T).tocoo()
    for start in range(0, len(gram.data), block):
        end = start+block
        gram.data[start:end] *= np.sum(codes[gram.row[start:end]]*codes[gram.col[start:end]], axis=1)
    result = gram.tocsr(); result.eliminate_zeros()
    return result


def debias_amplitudes(templates, codes, target, initial,
                      max_iterations=4000, tolerance=1e-6):
    """Solve fixed-template NNLS through a sparse Gram matrix, without a dense
    pixel-by-source design matrix. Signed prepared templates and spectral
    bases (e.g. GLS-whitened codes) are supported; amplitudes are nonnegative.
    """
    templates = csr_matrix(templates, dtype=float)
    codes = np.asarray(codes, float); target = np.asarray(target, float)
    initial = np.asarray(initial, float)
    n = templates.shape[0]
    if (codes.ndim != 2 or codes.shape[0] != n or target.shape != (templates.shape[1], codes.shape[1])
            or initial.shape != (n,) or np.any(initial < 0)
            or not all(np.isfinite(x).all() for x in [templates.data, codes, target, initial])
            or not isinstance(max_iterations, int) or max_iterations < 1
            or not np.isfinite(tolerance) or tolerance <= 0):
        raise ValueError('valid finite fixed-template amplitude problem required')
    if n == 0:
        return initial.copy(), dict(converged=True, iterations=0, fixed_point_relative=0.)
    gram = spectral_gram(templates, codes)
    b = np.sum((templates @ target)*codes, axis=1)
    lipschitz = float(np.max(np.asarray(abs(gram).sum(axis=1))))
    if lipschitz <= 0:
        return np.zeros(n), dict(converged=False, iterations=0, reason='unsupported_forward_templates')
    x = initial.copy(); z = x.copy(); t = 1.; relative = None
    for iteration in range(1, max_iterations+1):
        previous = x
        x = np.maximum(z-(gram @ z-b)/lipschitz, 0.)
        old_t = t; t = .5*(1+np.sqrt(1+4*t*t))
        z = x+(old_t-1)/t*(x-previous)
        if iteration % 25 == 0 or iteration == max_iterations:
            fixed = np.maximum(x-(gram @ x-b)/lipschitz, 0.)
            relative = float(np.linalg.norm(x-fixed)/max(np.linalg.norm(x), 1e-12))
            if relative <= tolerance:
                break
    return x, dict(converged=bool(relative <= tolerance), iterations=iteration,
                   fixed_point_relative=relative, lipschitz=lipschitz,
                   objective_without_constant=float(.5*(x @ (gram @ x))-b @ x),
                   nonnegative=True, penalty=0.)


def fit_debiased_fold_sources(stack, codebook, operator, heldout, penalty=.5,
                              iterations=150, minimum_norm=2.,
                              max_relative_error=.6, border=8):
    """Sparse selection followed by fitting-only amplitude refit and scoring."""
    state = fit_fold_sources(stack, codebook, operator, heldout, penalty=penalty,
                             iterations=iterations, minimum_norm=minimum_norm,
                             max_relative_error=max_relative_error, border=border)
    stack = np.asarray(stack, float); book = np.asarray(codebook, float)
    fitting = np.asarray(state['fitting_readouts']); held = np.asarray(state['heldout_readouts'])
    density = state['density']
    components = multiword_components(density, min_area=1, border=border)
    components = [c for c in components if c['peak']*np.linalg.norm(book[c['word']]) >= minimum_norm]
    if not components:
        state.update(source_rows=[], debias=dict(converged=True, iterations=0, sources=0))
        return state
    templates = []
    h, w = stack.shape[1:]
    for c in components:
        latent = np.zeros(h*w); latent[c['indices']] = c['coefficients']/c['peak']
        template = tiled_forward(latent.reshape(h, w, 1), operator)[:, :, 0]
        templates.append(csr_matrix(template.reshape(1, -1)))
    templates = vstack(templates, format='csr')
    labels = np.asarray([c['word'] for c in components])
    initial = np.asarray([c['peak'] for c in components])
    codes = book[labels]
    predicted_fit = tiled_forward(density @ book[:, fitting], operator).reshape(-1, len(fitting))
    nuisance_fit = predicted_fit-templates.T @ (initial[:, None]*codes[:, fitting])
    measured_fit = np.moveaxis(stack[fitting], 0, -1).reshape(-1, len(fitting))
    amplitudes, info = debias_amplitudes(templates, codes[:, fitting], measured_fit-nuisance_fit, initial)
    prediction = state['prediction'].reshape(-1, len(held))
    prediction = prediction+templates.T @ ((amplitudes-initial)[:, None]*codes[:, held])
    residual = np.moveaxis(stack[held], 0, -1).reshape(-1, len(held))-prediction
    projection = templates @ residual
    energies = np.asarray(templates.multiply(templates).sum(axis=1)).ravel()
    fit_norm = np.linalg.norm(book[:, fitting], axis=1)
    available = fit_norm > 0
    unit = book[available][:, fitting]/fit_norm[available, None]
    _, group, counts = np.unique(np.round(unit, 12), axis=0, return_inverse=True, return_counts=True)
    identified = np.zeros(len(book), bool); identified[available] = counts[group] == 1
    rows = []
    for i, c in enumerate(components):
        word = c['word']; expected = amplitudes[i]*codes[i, held]
        if not identified[word] or energies[i] <= 0 or np.linalg.norm(expected) <= 1e-12:
            continue
        measured = projection[i]/energies[i]+expected
        error = float(np.linalg.norm(measured-expected)/np.linalg.norm(expected))
        gain = float(energies[i]*(expected @ expected)+2*(projection[i] @ expected))
        rows.append(dict(word=word, position=c['position'], peak=float(amplitudes[i]),
                         penalized_peak=float(initial[i]),
                         heldout_on_readouts=int((book[word, held] > 0).sum()),
                         expected_trace=expected.tolist(), measured_trace=measured.tolist(),
                         relative_error=error, deletion_sse_gain=gain,
                         accepted=bool(info['converged'] and gain > 0 and error <= max_relative_error)))
    state.update(source_rows=rows, prediction=prediction.reshape(h, w, len(held)), debias=info,
                 interpretation='Sparse-selected fixed-source NNLS refit on fitting readouts only; withheld prediction, not molecular truth/FDR')
    return state
