"""Exact cross-label statistics for stochastically multiplexed localization data.

Each labelled molecule carries one of K labels (DNA-PAINT docking sequences, read
out in separate rounds), drawn independently with probabilities pi and kept for
the whole acquisition. Localizations inherit their molecule's label; molecule
identities are never used by the estimators below.

For a pair kernel phi, ordered localization pairs split into a same-molecule part
S and a between-molecule part W. Conditional on the configuration (positions and
counts) and averaging only over label draws,

    E[T_same]  = S + s2 W,        E[T_cross] = (1 - s2) W,

with s2 = sum(pi**2). Star triplets (a centre localization plus two neighbours,
all three labels distinct) involve three distinct molecules and satisfy

    E[T_star] = (1 - 3 s2 + 2 s3) W3,      s3 = sum(pi**3).
"""
import numpy as np
from scipy.spatial import cKDTree


def label_moments(pi):
    pi = np.asarray(pi, float)
    pi = pi / pi.sum()
    s2, s3 = float((pi ** 2).sum()), float((pi ** 3).sum())
    return dict(s2=s2, s3=s3, cross=1 - s2, distinct3=1 - 3 * s2 + 2 * s3)


def pair_counts(xy, label, edges, K):
    """Ordered localization pairs (i != j) by label pair and distance bin.

    Returns an array N[a, b, bin] with N[a, b] == N[b, a].
    """
    edges = np.asarray(edges, float)
    trees = [cKDTree(xy[label == a]) for a in range(K)]
    sizes = [int((label == a).sum()) for a in range(K)]
    N = np.zeros((K, K, len(edges) - 1))
    for a in range(K):
        for b in range(a, K):
            if sizes[a] == 0 or sizes[b] == 0:
                continue
            cum = trees[a].count_neighbors(trees[b], edges).astype(float)
            if a == b:
                cum -= sizes[a]  # remove i == j
            h = np.diff(cum)
            N[a, b] = h
            N[b, a] = h
    return N


def split_pairs(N):
    K = N.shape[0]
    same = sum(N[a, a] for a in range(K))
    cross = N.sum((0, 1)) - same
    return same, cross


def decompose(N, pi):
    """Unbiased same-molecule (S) and between-molecule (W) pair histograms."""
    m = label_moments(pi)
    same, cross = split_pairs(N)
    W = cross / m['cross']
    return dict(S=same - m['s2'] * W, W=W, same=same, cross=cross)


def star_counts(xy, label, radius, K):
    """Per-localization neighbour counts by label within `radius` (self excluded)."""
    C = np.zeros((len(xy), K))
    for b in range(K):
        pts = xy[label == b]
        if len(pts):
            C[:, b] = cKDTree(pts).query_ball_point(xy, radius, return_length=True)
    C[np.arange(len(xy)), label] -= 1  # each point finds itself in its own label
    return C


def star_triplets(xy, label, radius, K, counts=None):
    """Sum over centres i of ordered neighbour pairs (j, k) with three distinct labels."""
    C = star_counts(xy, label, radius, K) if counts is None else counts
    C = C.copy()
    C[np.arange(len(xy)), label] = 0
    return float((C.sum(1) ** 2 - (C ** 2).sum(1)).sum())


def lens_area(r, R):
    """Overlap area of two discs of radius R whose centres are r apart."""
    r = np.minimum(np.asarray(r, float), 2 * R)
    return 2 * R ** 2 * np.arccos(r / (2 * R)) - r / 2 * np.sqrt(4 * R ** 2 - r ** 2)


def lens_integral(edges, h, R):
    """Integral over two discs of radius R of h(|r1 - r2|), from binned h:
    sum over bins of h * lens_area * bin area (needs edges[-1] >= 2R)."""
    edges = np.asarray(edges, float)
    c = (edges[:-1] + edges[1:]) / 2
    ann = np.pi * np.diff(edges ** 2)
    inside = edges[1:] <= 2 * R + 1e-9
    return (np.asarray(h)[..., inside] * lens_area(c[inside], R) * ann[inside]).sum(-1)


def connected_triplets(xy, label, K, radius, annulus, L):
    """Per-centre connected triplet terms with three distinct labels.

    For centre i (label a), e_b = n_b - lambda_b pi R^2 is the excess count of
    label-b points within `radius`, with lambda_b the label-b density in
    `annulus` around i. Returns q_i = sum over ordered label pairs (b, c),
    b != c, both != a, of e_b e_c - lambda_b lambda_c L[b, c], where L[b, c] is
    the two-disc integral of the between-molecule pair excess h_bc.

    For complexes placed at random with independent labels, the expectation of
    q_i is twice the probability that i has two labelled complex partners within
    `radius` carrying labels distinct from each other and from a. It is zero for
    any monomer/dimer mixture, at any labelling efficiency and density.
    """
    A_in = np.pi * radius ** 2
    A_ann = np.pi * (annulus[1] ** 2 - annulus[0] ** 2)
    n_in = star_counts(xy, label, radius, K)
    n_ann = star_counts(xy, label, annulus[1], K) - star_counts(xy, label, annulus[0], K)
    lam = n_ann / A_ann
    e = n_in - lam * A_in
    e[np.arange(len(xy)), label] = 0
    lam_off = lam.copy()
    lam_off[np.arange(len(xy)), label] = 0
    # sum_{b != c} e_b e_c = (sum e)^2 - sum e^2 over labels other than a
    prod = e.sum(1) ** 2 - (e ** 2).sum(1)
    pair = np.einsum('ib,bc,ic->i', lam_off, L, lam_off) - np.einsum('ib,bb,ib->i', lam_off, L, lam_off)
    return prod - pair


def lens_angle_profile(edges, h, R, angle_edges, rng, n=400000):
    """Angle-resolved two-disc integral of h(|r1 - r2|): for r1, r2 uniform in a
    disc of radius R around the origin, A^2 E[h(|r1 - r2|) 1(angle(r1, r2) in bin)].
    `h` has a leading label-pair shape (..., bins); returns (..., angle bins)."""
    edges = np.asarray(edges, float)
    A = np.pi * R ** 2
    rad = R * np.sqrt(rng.random((2, n)))
    phi = rng.uniform(0, 2 * np.pi, (2, n))
    p = np.stack([rad * np.cos(phi), rad * np.sin(phi)], -1)          # (2, n, 2)
    d = np.hypot(*(p[0] - p[1]).T)
    ang = np.abs((phi[0] - phi[1] + np.pi) % (2 * np.pi) - np.pi)       # in [0, pi]
    db = np.clip(np.digitize(d, edges) - 1, 0, len(edges) - 2)
    ab = np.clip(np.digitize(ang, angle_edges) - 1, 0, len(angle_edges) - 2)
    h = np.asarray(h)
    flat = h.reshape(-1, h.shape[-1])
    out = np.zeros((flat.shape[0], len(angle_edges) - 1))
    for k, row in enumerate(flat):
        out[k] = np.bincount(ab, weights=row[db], minlength=len(angle_edges) - 1) / n * A ** 2
    return out.reshape(h.shape[:-1] + (len(angle_edges) - 1,))


def connected_triplet_angles(xy, label, K, radius, annulus, L_angle, angle_edges):
    """Angle-resolved connected triplets: for each centre i (label a) and each
    ordered pair of neighbours (j, k) within `radius` with labels b != c, both != a,
    the angle j-i-k is binned. Background terms (uniform in angle) and the
    neighbour-neighbour pair term L_angle[b, c, bin] are subtracted as in
    `connected_triplets`; summing the result over angle bins gives that statistic.
    Returns the per-centre (n, bins) array."""
    nb_bins = len(angle_edges) - 1
    width = np.diff(angle_edges) / np.pi
    A_in = np.pi * radius ** 2
    A_ann = np.pi * (annulus[1] ** 2 - annulus[0] ** 2)
    n_in = star_counts(xy, label, radius, K)
    lam = (star_counts(xy, label, annulus[1], K) - star_counts(xy, label, annulus[0], K)) / A_ann
    tree = cKDTree(xy)
    nbrs = tree.query_ball_point(xy, radius)
    out = np.zeros((len(xy), nb_bins))
    for i, js in enumerate(nbrs):
        a = label[i]
        js = np.asarray([j for j in js if j != i and label[j] != a], int)
        if len(js) >= 2:
            v = xy[js] - xy[i]
            th = np.arctan2(v[:, 1], v[:, 0])
            lj = label[js]
            jj, kk = np.meshgrid(np.arange(len(js)), np.arange(len(js)), indexing='ij')
            ok = lj[jj] != lj[kk]
            ang = np.abs((th[jj] - th[kk] + np.pi) % (2 * np.pi) - np.pi)[ok]
            out[i] = np.histogram(ang, angle_edges)[0]
        # background expectation, uniform in angle: sum_{b != c} (bg_c n_b + bg_b n_c - bg_b bg_c)
        bg = lam[i] * A_in
        bg[a] = 0
        nn = n_in[i].copy(); nn[a] = 0
        tot = 2 * (bg.sum() * nn.sum() - (bg * nn).sum()) - (bg.sum() ** 2 - (bg ** 2).sum())
        lo = lam[i].copy(); lo[a] = 0
        pair = np.einsum('b,bcm,c->m', lo, L_angle, lo) - np.einsum('b,bbm,b->m', lo, L_angle, lo)
        out[i] -= tot * width + pair
    return out


# ---------------------------------------------------------------------------
# Oracle quantities (molecule identities known); used only to verify the exact
# identities in simulation.

def oracle_pairs(xy, mol, edges):
    """Ordered same-molecule (S) and between-molecule (W) pair histograms, plus
    per-molecule-pair counts within each bin for the exact label variance."""
    edges = np.asarray(edges, float)
    tree = cKDTree(xy)
    i, j = tree.query_pairs(edges[-1], output_type='ndarray').T
    d = np.hypot(*(xy[i] - xy[j]).T)
    b = np.digitize(d, edges) - 1
    ok = (b >= 0) & (b < len(edges) - 1)
    i, j, b = i[ok], j[ok], b[ok]
    same = mol[i] == mol[j]
    nb = len(edges) - 1
    S = 2 * np.bincount(b[same], minlength=nb).astype(float)
    W = 2 * np.bincount(b[~same], minlength=nb).astype(float)
    lo, hi = np.minimum(mol[i], mol[j])[~same], np.maximum(mol[i], mol[j])[~same]
    return S, W, (lo, hi, b[~same])


def exact_cross_variance(pairs, pi, nbins, nmol):
    """Exact label-draw variance of T_cross in each bin, given the configuration.

    T_cross = sum_{s<t} v_st (1 - 1[c_s = c_t]) with v_st the ordered pair count.
    Indicators of disjoint molecule pairs are independent; pairs sharing one
    molecule have covariance s3 - s2**2 (zero for uniform labels).
    """
    m = label_moments(pi)
    var1, cov1 = m['s2'] - m['s2'] ** 2, m['s3'] - m['s2'] ** 2
    lo, hi, b = pairs
    out = np.zeros(nbins)
    for k in range(nbins):
        sel = b == k
        if not sel.any():
            continue
        key = lo[sel] * nmol + hi[sel]
        u, v = np.unique(key, return_counts=True)
        v = 2.0 * v  # ordered pairs
        s, t = u // nmol, u % nmol
        total = var1 * (v ** 2).sum()
        if cov1 != 0:
            # sum over molecules of (sum_t v_st)^2 - sum_t v_st^2, i.e. pairs sharing s
            per = np.bincount(np.r_[s, t], weights=np.r_[v, v], minlength=nmol)
            per2 = np.bincount(np.r_[s, t], weights=np.r_[v, v] ** 2, minlength=nmol)
            total += cov1 * (per ** 2 - per2).sum()
        out[k] = total
    return out


def oracle_star(xy, mol, radius):
    """Star-triplet sum over three distinct molecules (labels ignored)."""
    tree = cKDTree(xy)
    nb = tree.query_ball_point(xy, radius)
    total = 0.0
    for i, js in enumerate(nb):
        js = np.asarray(js)
        js = js[mol[js] != mol[i]]
        if len(js) < 2:
            continue
        _, c = np.unique(mol[js], return_counts=True)
        total += c.sum() ** 2 - (c ** 2).sum()
    return total
