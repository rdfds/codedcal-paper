"""Mixture-preserving components and checks against a fitted forward model.

This is an in-sample model check, not independent molecular identity evidence.
Other fitted sources are removed before projecting a candidate's measured trace.
Exact same-template code aliases cause abstention; arbitrary mixture uniqueness
is not certified. The old winner-only/native-equivalence adapter is unchanged.
"""
import numpy as np
from scipy.ndimage import label
from scipy.optimize import nnls
from scipy.sparse import csr_matrix, vstack


def multiword_components(density, threshold=0., min_area=1, border=8):
    """Keep connected positive components separately for every word.

    Positions use coefficient-weighted camera-pixel centroids. Support indices
    are C-order flat indices in the spatial grid; no upsampling is assumed.
    """
    density = np.asarray(density, float)
    if (density.ndim != 3 or not np.isfinite(density).all()
            or np.any(density < 0) or not np.isfinite(threshold) or threshold < 0
            or not isinstance(min_area, (int, np.integer)) or min_area < 1
            or not np.isfinite(border) or border < 0):
        raise ValueError('finite nonnegative density and valid component settings required')
    h, w, _ = density.shape
    components = []
    for word in range(density.shape[2]):
        image = density[:, :, word]
        ids, _ = label(image > threshold)
        sizes = np.bincount(ids.ravel())
        for component in np.flatnonzero(sizes >= min_area):
            if component == 0:
                continue
            yy, xx = np.where(ids == component)
            weights = image[yy, xx]
            mass = weights.sum()
            x, y = float(weights @ xx / mass), float(weights @ yy / mass)
            if border <= x < w-border and border <= y < h-border:
                components.append(dict(word=word, position=[x, y], area=len(xx),
                                       indices=np.ravel_multi_index((yy, xx), (h, w)),
                                       coefficients=weights, peak=float(weights.max())))
    return components


def tiled_forward(coefficients, operator):
    """Apply the exact pixel-grid tile operator used in sparse fitting.

    The caller must pass the same preparation/PSF operator and boundary model.
    Last-axis columns can be words, readouts, or one component's spatial shape.
    """
    x = np.asarray(coefficients, float)
    a = np.asarray(operator, float)
    side = int(np.sqrt(a.shape[0])) if a.ndim == 2 else 0
    if (x.ndim != 3 or side < 1 or a.shape != (side*side, side*side)
            or x.shape[0] % side or x.shape[1] % side
            or not np.isfinite(x).all() or not np.isfinite(a).all()):
        raise ValueError('finite pixel-grid square tiles required')
    result = np.zeros_like(x)
    for y in range(0, x.shape[0], side):
        for z in range(0, x.shape[1], side):
            tile = x[y:y+side, z:z+side]
            if not np.any(tile):
                continue
            result[y:y+side, z:z+side] = (a @ tile.reshape(side*side, -1, order='F')).reshape(tile.shape, order='F')
    return result


def joint_source_consistency(stack, density, codebook, operator,
                             threshold=0., min_area=1, border=8,
                             max_distance=.6, minimum_norm=2.,
                             alias_tolerance=1e-8):
    """Check decontaminated traces and component-deletion prediction evidence.

    Labels and coefficients are never re-estimated. Components with nonpositive
    deletion evidence are pruned before the photometric check, so an invented
    source cannot contaminate its neighbor's check. An accepted source must
    have positive deletion evidence, nonzero measured contrast, and agree with
    its word after subtracting other fitted contributions. These are model
    diagnostics on the fitting data; positivity is not a statistical discovery
    threshold. Minimum norm applies to both the fitted and measured source
    vectors; it screens negligible fitted components before allocating templates.
    It uses forward-model coefficient units, not the old
    rectified-photon units. It needs task-specific development qualification.
    """
    y = np.asarray(stack, float)
    density = np.asarray(density, float)
    book = np.asarray(codebook, float)
    if (y.ndim != 3 or density.ndim != 3 or book.ndim != 2
            or density.shape[:2] != y.shape[1:]
            or book.shape != (density.shape[2], y.shape[0])
            or not np.isfinite(y).all() or not np.isfinite(book).all()
            or np.any(book < 0) or np.any(np.linalg.norm(book, axis=1) == 0)):
        raise ValueError('finite matching stack/density/codebook required')
    if (not np.isfinite(max_distance) or not 0 <= max_distance <= 2
            or not np.isfinite(minimum_norm) or minimum_norm < 0
            or not np.isfinite(alias_tolerance) or not 0 < alias_tolerance < 1):
        raise ValueError('invalid quality settings')
    components = multiword_components(density, threshold, min_area, border)
    original_count = len(components)
    components = [c for c in components
                  if c['peak']*np.linalg.norm(book[c['word']]) >= minimum_norm]
    prediction = tiled_forward(density @ book, operator)
    residual = np.moveaxis(y, 0, -1) - prediction
    templates = []
    h, w = density.shape[:2]
    for c in components:
        latent = np.zeros(h*w)
        latent[c['indices']] = c['coefficients'] / c['peak']
        template = tiled_forward(latent.reshape(h, w, 1), operator)[:, :, 0]
        templates.append(csr_matrix(template.reshape(1, -1)))
    n = len(components)
    distances = np.zeros(n); norms = np.zeros(n); evidence = np.zeros(n)
    traces = np.zeros((n, len(y))); pruned = np.zeros(n, bool)
    ambiguous = np.zeros(n, bool)
    unit_book = book / np.linalg.norm(book, axis=1)[:, None]
    # Restrict the exact-alias test to effectively identical spatial templates.
    # Other positional/near-collinear aliases are an explicit unresolved scope.
    if templates:
        templates = vstack(templates, format='csr')
        energies = np.sqrt(np.asarray(templates.multiply(templates).sum(axis=1)).ravel())
        assigned = book[[c['word'] for c in components]]
        strengths = np.array([c['peak'] for c in components])[:, None] * assigned
        flat_residual = residual.reshape(-1, len(y))
        for _ in range(n+1):
            projection = templates @ flat_residual
            gains = energies**2*np.sum(strengths**2, axis=1) + 2*np.sum(projection*strengths, axis=1)
            evidence[~pruned] = gains[~pruned]
            remove = ~pruned & ((gains <= 0) | (energies == 0))
            if not remove.any():
                break
            flat_residual += templates[remove].T @ strengths[remove]
            pruned |= remove
        traces = projection / np.maximum(energies[:, None]**2, np.finfo(float).tiny)
        traces += strengths * ~pruned[:, None]
        norms = np.linalg.norm(traces, axis=1)
        distances = np.linalg.norm(traces / np.maximum(norms[:, None], np.finfo(float).tiny)
                                   - assigned / np.linalg.norm(assigned, axis=1)[:, None], axis=1)
        unit_templates = templates.multiply(1/np.maximum(energies[:, None], np.finfo(float).tiny))
        similarities = (unit_templates @ unit_templates.T).tocsr()
        visited = set()
        for i, c in enumerate(components):
            if pruned[i]:
                continue
            row = similarities.getrow(i)
            group = row.indices[row.data >= 1-alias_tolerance]
            group = group[~pruned[group]]
            key = tuple(group)
            if key in visited:
                continue
            visited.add(key)
            # Scale each source by its template norm for a common spatial shape.
            target = sum(components[j]['peak'] * energies[j] * book[components[j]['word']]
                         for j in group)
            words = {components[j]['word'] for j in group}
            others = np.array([k for k in range(len(book)) if k not in words], int)
            # If positive active amplitudes already lie in a dependent code
            # set, a small null-space change preserves positivity and the
            # observation. Excluding all active words from the alternative
            # search must not conceal this within-set ambiguity.
            active_codes = book[[components[j]['word'] for j in group]]
            alias = bool(len(group) > 1 and
                         np.linalg.matrix_rank(active_codes, tol=alias_tolerance) < len(group))
            if len(group) == 1:
                alias = bool(len(others) and np.any(np.linalg.norm(unit_book[others]-unit_book[c['word']], axis=1) <= alias_tolerance))
            elif not alias and len(others) and np.linalg.norm(target) > 0:
                try:
                    _, error = nnls(book[others].T, target, maxiter=max(100, 20*len(others)))
                    alias = error <= alias_tolerance*np.linalg.norm(target)
                except RuntimeError:
                    # An unresolved test is unsupported, never a positive call.
                    alias = True
            ambiguous[group] |= alias
    keep = (norms > 0) & (norms >= minimum_norm) & (distances <= max_distance) & (evidence > 0) & ~ambiguous & ~pruned
    return components, keep, dict(distance=distances.tolist(), norm=norms.tolist(),
                                 deletion_sse_gain=evidence.tolist(), ambiguous=ambiguous.tolist(),
                                 pruned_nonpositive_evidence=pruned.tolist(),
                                 measured_decontaminated_trace=traces.tolist(),
                                 original_components=original_count,
                                 screened_low_fitted_norm=original_count-len(components),
                                 residual_sse=float(np.sum(residual**2)),
                                 interpretation='In-sample forward-model check; not molecular truth, FDR or arbitrary-mixture identifiability')
