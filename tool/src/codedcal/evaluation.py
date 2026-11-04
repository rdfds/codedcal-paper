"""Identity and composition evaluation with explicit reference limitations.

Reference labels are truth only when the caller has independently established
them. A released decode or a repeated observation supplies a consistency
reference. None of these functions estimates molecular FDR from blank counts.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from scipy.spatial.distance import jensenshannon
from scipy.ndimage import affine_transform


def validate_calls(xy, words, n_words):
    xy, words = np.asarray(xy, dtype=float), np.asarray(words)
    if xy.shape != (len(words), 2) or words.ndim != 1:
        raise ValueError('calls must have matching (n,2) coordinates and word IDs')
    if not np.isfinite(xy).all() or not np.issubdtype(words.dtype, np.integer):
        raise ValueError('call coordinates must be finite; word IDs integer')
    if n_words < 2 or np.any((words < 0) | (words >= n_words)):
        raise ValueError('word ID outside codebook')
    return xy, words


def spatial_match(reference, prediction, radius=1.5, max_component=2000):
    """Label-independent maximum-cardinality, then minimum-distance matching.

    Resolve only connected spatial components. Components exceeding the memory
    limit remain unmatched and are disclosed; they never disappear from scoring.
    Report ambiguity rather than treating a spatial match as physical identity.
    """
    a, b = np.asarray(reference, float), np.asarray(prediction, float)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != 2 or b.shape[1] != 2:
        raise ValueError('coordinates must be (n,2)')
    if not np.isfinite(a).all() or not np.isfinite(b).all() or not np.isfinite(radius) or radius <= 0:
        raise ValueError('coordinates and positive radius must be finite')
    if not isinstance(max_component, int) or max_component < 2:
        raise ValueError('max_component must be integer >=2')
    if not len(a) or not len(b):
        return np.empty((0, 2), int), dict(ambiguous_reference=0, ambiguous_prediction=0,
                                         skipped_components=0, candidate_edges=0)
    edges = cKDTree(a).query_ball_tree(cKDTree(b), radius)
    parent = np.arange(len(a) + len(b))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return int(i)
    degree_b = np.zeros(len(b), int)
    for i, neighbors in enumerate(edges):
        for j in neighbors:
            parent[find(len(a) + j)] = find(i)
            degree_b[j] += 1
    groups = {}
    for i, neighbors in enumerate(edges):
        if neighbors:
            groups.setdefault(find(i), [[], set()])[0].append(i)
            groups[find(i)][1].update(neighbors)
    pairs, skipped = [], 0
    for aa, bb in groups.values():
        bb = sorted(bb)
        nr, nc = len(aa), len(bb)
        if nr + nc > max_component:
            skipped += 1
            continue
        dist = np.linalg.norm(a[aa, None, :] - b[None, bb, :], axis=2)
        reward = (min(nr, nc) + 1) * radius
        cost = np.zeros((nr + nc, nr + nc))
        cost[:nr, :nc] = np.where(dist <= radius, dist - reward, reward)
        rr, cc = linear_sum_assignment(cost)
        for i, j in zip(rr, cc):
            if i < nr and j < nc and dist[i, j] <= radius:
                pairs.append((aa[i], bb[j]))
    return np.asarray(sorted(pairs), dtype=int).reshape(-1, 2), dict(
        ambiguous_reference=sum(len(x) > 1 for x in edges),
        ambiguous_prediction=int((degree_b > 1).sum()), skipped_components=skipped,
        candidate_edges=sum(map(len, edges)))


def evaluate_calls(reference_xy, reference_words, call_xy, call_words,
                   is_blank, radius=1.5, reference_kind='simulation_truth'):
    """Score all calls; distinguish identity errors from composition errors.

    Identity F1 counts blanks as erroneous calls. Gene precision additionally
    reports the corresponding gene-only denominator. Empty risks are undefined.
    """
    blank = np.asarray(is_blank)
    if blank.ndim != 1 or blank.dtype != bool or not np.any(blank) or not np.any(~blank):
        raise ValueError('is_blank must be a boolean book mask with both classes')
    if reference_kind not in ('simulation_truth', 'unperturbed_call_reference', 'repeat_consistency'):
        raise ValueError('reference_kind must explicitly identify what is scored')
    a, aw = validate_calls(reference_xy, reference_words, len(blank))
    b, bw = validate_calls(call_xy, call_words, len(blank))
    if np.any(blank[aw]):
        raise ValueError('gene reference cannot contain blank words')
    pairs, diagnostic = spatial_match(a, b, radius)
    correct_mask = aw[pairs[:, 0]] == bw[pairs[:, 1]]
    correct = int(correct_mask.sum())
    counts = np.bincount(bw, minlength=len(blank))
    refcounts = np.bincount(aw, minlength=len(blank))
    gene, blanks = int(counts[~blank].sum()), int(counts[blank].sum())
    r, c = len(a), len(b)
    tv = js = None
    if r and gene:
        p, q = refcounts[~blank] / r, counts[~blank] / gene
        tv = float(np.abs(p-q).sum()/2)
        js = float(jensenshannon(p, q, base=2)**2)
    distance = np.linalg.norm(a[pairs[:, 0]] - b[pairs[:, 1]], axis=1)
    return dict(reference_kind=reference_kind, n_reference=r, calls=c, gene_calls=gene,
                blank_calls=blanks, spatial_matches=len(pairs), correct_gene_calls=correct,
                wrong_word_matches=int(len(pairs)-correct), unmatched_calls=c-len(pairs),
                missed_reference=r-len(pairs), identity_precision=correct/c if c else None,
                gene_precision=correct/gene if gene else None,
                identity_recall=correct/r if r else None,
                identity_f1=2*correct/(r+c) if r+c else None,
                false_call_fraction=(c-correct)/c if c and reference_kind == 'simulation_truth' else None,
                blank_proxy=blanks/gene*int((~blank).sum())/int(blank.sum()) if gene else None,
                composition_tv=tv, composition_js=js,
                matched_rms_px=float(np.sqrt(np.mean(distance**2))) if len(distance) else None,
                matching=diagnostic)


def identity_affines(n_bits):
    return np.repeat(np.array([[1., 0.], [0., 1.], [0., 0.]])[None], n_bits, axis=0)


def compose_affines(first, second):
    """Maps use row coordinates: raw @ map[:2] + map[2]."""
    a, b = np.asarray(first, float), np.asarray(second, float)
    if a.shape != b.shape or a.ndim != 3 or a.shape[1:] != (3, 2) or not np.isfinite([a,b]).all():
        raise ValueError('maps must be matching finite (bits,3,2) arrays')
    return np.concatenate([a[:, :2] @ b[:, :2],
                           ((a[:, 2, None, :] @ b[:, :2])[:, 0] + b[:, 2])[:, None, :]], axis=1)


def warp_readouts(images, raw_to_aligned, origin=(0., 0.)):
    """Shared inverse resampler; interpolation loss is not measurement truth."""
    images, maps = np.asarray(images), np.asarray(raw_to_aligned, float)
    origin = np.asarray(origin, float)
    if images.ndim != 3 or maps.shape != (len(images), 3, 2) or origin.shape != (2,):
        raise ValueError('images (bits,y,x), maps (bits,3,2), origin (x,y) required')
    if not np.isfinite(images).all() or not np.isfinite(maps).all() or not np.isfinite(origin).all():
        raise ValueError('images, maps and origin must be finite')
    out = np.empty(images.shape, np.float32)
    for k, (image, coef) in enumerate(zip(images, maps)):
        if abs(np.linalg.det(coef[:2])) < 1e-9:
            raise ValueError('singular affine map')
        inv = np.linalg.inv(coef[:2])
        offset = (origin-coef[2]) @ inv-origin
        out[k] = affine_transform(image.astype(np.float32), inv.T[::-1, ::-1],
                                  offset[::-1], order=1, mode='constant', cval=0,
                                  prefilter=False)
    return out
