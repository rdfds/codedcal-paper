"""Cross-label correlation of point sets with random permanent labels.

Points are localizations or molecule positions (nm) with integer labels
0..n_labels-1. Pair counts are accumulated per spatial block so that
uncertainty can be estimated by resampling blocks. The cross-label pair
correlation g(r) is normalized by a null in which labels are displaced from one
another by micrometres, corrected on a plateau annulus for the overlap the
displacement loses.
"""
import numpy as np
from scipy.spatial import cKDTree

from .labels import (connected_triplet_angles, connected_triplets, label_moments, lens_angle_profile,
                     lens_integral)

EDGES = np.r_[np.arange(0, 60, 1.0), np.arange(60, 101, 5.0)]   # nm
PLATEAU = np.array([300.0, 600.0])       # nm, normalization annulus
PLATEAU_FRACTION = 0.05                  # source subsample for plateau counts
NULL_SMOOTH = (20.0, 100.0)              # nm, annulus defining the null pair density
BLOCK = 2000.0                           # nm, spatial resampling unit
BASELINE = (60.0, 100.0)                 # nm, baseline window for excess statistics
TRIPLET_RADIUS = 30.0                    # nm
TRIPLET_ANNULUS = (100.0, 250.0)         # nm, local background for triplets
ANGLE_EDGES = np.linspace(0, np.pi, 13)  # 15-degree bins


def _n_labels(label, n_labels):
    return int(label.max()) + 1 if n_labels is None else int(n_labels)


def blocks(xy, size=BLOCK):
    """Index of the size x size nm block containing each point."""
    ij = np.floor((xy - xy.min(0)) / size).astype(int)
    _, idx = np.unique(ij, axis=0, return_inverse=True)
    return idx.ravel()


def block_pair_counts(xy, label, edges=EDGES, block=None, shift=None, n_labels=None):
    """Ordered pair counts N[k, a, b, bin] with the first point in block k, and
    point counts n[k, a].

    With `shift`, every label's points are displaced by shift[label] (null
    construction: relative displacements shift[b] - shift[a]); the self-pair
    correction is then not applied.
    """
    K = _n_labels(label, n_labels)
    block = blocks(xy) if block is None else block
    nb = block.max() + 1
    targets = []
    for b in range(K):
        pts = xy[label == b]
        if shift is not None:
            pts = pts + shift[b]
        targets.append(cKDTree(pts))
    N = np.zeros((nb, K, K, len(edges) - 1))
    n = np.zeros((nb, K))
    for k in range(nb):
        for a in range(K):
            sel = (block == k) & (label == a)
            n[k, a] = sel.sum()
            if not sel.any():
                continue
            src = cKDTree(xy[sel] if shift is None else xy[sel] + shift[a])
            for b in range(K):
                cum = src.count_neighbors(targets[b], edges).astype(float)
                if a == b and shift is None:
                    cum -= sel.sum()
                N[k, a, b] = np.diff(cum)
    return N, n


def cross_g(N, n):
    """Per-label-pair and pooled cross-label pair counts normalized by sum n_a n_b."""
    K = N.shape[1]
    Ns, ns = N.sum(0), n.sum(0)
    G = np.zeros_like(Ns)
    for a in range(K):
        for b in range(K):
            G[a, b] = Ns[a, b] / max(ns[a] * ns[b], 1)
    off = ~np.eye(K, dtype=bool)
    pooled = Ns[off].sum(0) / (np.outer(ns, ns)[off].sum())
    return G, pooled


def decomposition(N, n):
    """Same-molecule pair histogram per label: S_a = N_aa - n_a^2 * pooled cross."""
    K = N.shape[1]
    Ns, ns = N.sum(0), n.sum(0)
    _, pooled = cross_g(N, n)
    S = np.array([Ns[a, a] - ns[a] ** 2 * pooled for a in range(K)])
    return S, pooled


def cross_annuli(xy, label, seed, shift=None, fraction=PLATEAU_FRACTION,
                 annuli=(NULL_SMOOTH, tuple(PLATEAU)), n_labels=None):
    """Pooled cross-label pair counts per annulus from a fixed source subsample,
    normalized by sum_{a != b} n_a(subsample) n_b like `cross_g`."""
    K = _n_labels(label, n_labels)
    rng = np.random.default_rng(seed)
    radii = np.array([r for a in annuli for r in a], float)
    counts = np.zeros(len(annuli)); norm = 0.0
    tgts = [xy[label == b] if shift is None else xy[label == b] + shift[b] for b in range(K)]
    tgt_trees = [cKDTree(t) for t in tgts]
    for a in range(K):
        src = xy[label == a]
        src = src[rng.random(len(src)) < fraction]
        if shift is not None:
            src = src + shift[a]
        ts = cKDTree(src)
        for b in range(K):
            if b == a:
                continue
            cum = ts.count_neighbors(tgt_trees[b], radii).astype(float)
            counts += cum[1::2] - cum[0::2]
            norm += len(src) * len(tgts[b])
    return counts / norm


def null_shifts(rng, n_shift=6, lo=1500.0, hi=3000.0, n_labels=4):
    """Labels displaced to evenly spaced directions at a common random radius, so
    every pair of labels is separated by at least 2 lo sin(pi / n_labels)."""
    out = []
    for _ in range(n_shift):
        ang = rng.uniform(0, 2 * np.pi) + 2 * np.pi * np.arange(n_labels) / n_labels
        out.append(rng.uniform(lo, hi) * np.column_stack([np.cos(ang), np.sin(ang)]))
    return out


def null_normalization(xy, label, rng, n_shift=6, n_labels=None):
    """Null pair density and plateau ratio from shifted configurations.

    Returns (null_density, plateau_ratio): the shifted-null cross pair density per
    unit annulus area over NULL_SMOOTH (pairs / sum n_a n_b / nm^2), and the
    observed/null ratio of pair counts on the plateau annulus.
    """
    K = _n_labels(label, n_labels)
    seed = int(rng.integers(1 << 31))
    obs = cross_annuli(xy, label, seed, n_labels=K)
    null = np.mean([cross_annuli(xy, label, seed, shift=s, n_labels=K)
                    for s in null_shifts(rng, n_shift, n_labels=K)], 0)
    area = np.pi * (NULL_SMOOTH[1] ** 2 - NULL_SMOOTH[0] ** 2)
    return null[0] / area, obs[1] / null[1]


def g_of_r(pooled, null_density, plateau_ratio, edges=EDGES):
    """Cross-label g(r) = observed pair density / (null density x plateau ratio).

    At r <= 100 nm edge losses differ from the null by under 1%, so the null
    enters through its mean density over NULL_SMOOTH; dividing by sparse null
    counts at small r would bias the ratio upward.
    """
    ann = np.pi * np.diff(edges ** 2)
    return pooled / (ann * null_density) / plateau_ratio


def point_pair_g(xy, label, rng, n_shift=12, fraction=1.0, block_size=BLOCK, n_labels=None):
    """Cross-label correlation for one point per molecule (or any modest point set):
    block counts, shifted-label normalization, pooled and per-label-pair g."""
    K = _n_labels(label, n_labels)
    block = blocks(xy, block_size)
    N, n = block_pair_counts(xy, label, block=block, n_labels=K)
    seed = int(rng.integers(1 << 31))
    obs = cross_annuli(xy, label, seed, fraction=fraction, n_labels=K)
    null = np.mean([cross_annuli(xy, label, seed, shift=s, fraction=fraction, n_labels=K)
                    for s in null_shifts(rng, n_shift, n_labels=K)], 0)
    nd = null[0] / (np.pi * (NULL_SMOOTH[1] ** 2 - NULL_SMOOTH[0] ** 2))
    pr = obs[1] / null[1]
    G, pooled = cross_g(N, n)
    pair_g = np.array([[g_of_r(G[a, b], nd, pr) for b in range(K)] for a in range(K)])
    return dict(block=block, N=N, n=n, null_density=nd, plateau_ratio=pr,
                g=g_of_r(pooled, nd, pr), pair_g=pair_g)


def baseline_excess(g, edges=EDGES, base=BASELINE):
    """g minus its area-weighted mean over `base` (last axis = distance bins)."""
    c = (edges[:-1] + edges[1:]) / 2
    ann = np.pi * np.diff(edges ** 2)
    sel = (c >= base[0]) & (c < base[1])
    b = (g[..., sel] * ann[sel]).sum(-1) / ann[sel].sum()
    return g - b[..., None]


def bootstrap(fn, N, n, rng, B=200):
    """Resample spatial blocks; fn(N, n) -> array."""
    nb = N.shape[0]
    out = []
    for _ in range(B):
        k = rng.integers(0, nb, nb)
        out.append(fn(N[k], n[k]))
    return np.array(out)


def registration_offsets(xy, label, radius=25.0, fraction=0.05, seed=0, n_labels=None):
    """Mean cross-label displacement (b minus a) among close pairs, per label pair,
    from a source subsample, with approximate standard errors and pair counts.

    Between-molecule displacements are symmetric in expectation, so a nonzero mean
    measures a relative registration offset between labels (imaging rounds).
    """
    K = _n_labels(label, n_labels)
    rng = np.random.default_rng(seed)
    out = np.zeros((K, K, 2)); se = np.zeros((K, K, 2)); cnt = np.zeros((K, K))
    pts = [xy[label == b] for b in range(K)]
    trees = [cKDTree(p) for p in pts]
    for a in range(K):
        src = pts[a][rng.random(len(pts[a])) < fraction]
        for b in range(a + 1, K):
            nb = trees[b].query_ball_point(src, radius)
            ii = np.repeat(np.arange(len(nb)), [len(x) for x in nb])
            if not len(ii):
                continue
            jj = np.concatenate([np.asarray(x, int) for x in nb])
            dvec = pts[b][jj] - src[ii]
            out[a, b] = dvec.mean(0); out[b, a] = -out[a, b]
            # approximate s.e. treating per-source mean displacements as independent
            per = np.array([np.bincount(ii, weights=dvec[:, c], minlength=len(src)) for c in range(2)]).T
            m = np.bincount(ii, minlength=len(src))
            ok = m > 0
            se[a, b] = se[b, a] = (per[ok] / m[ok, None]).std(0) / np.sqrt(ok.sum())
            cnt[a, b] = cnt[b, a] = len(ii)
    return out, se, cnt


def label_translations(offsets):
    """Least-squares per-label translations t (sum t = 0) from pairwise offsets
    o_ab = t_b - t_a. Returns translations and residual rms."""
    K = offsets.shape[0]
    rows, rhs = [], []
    for a in range(K):
        for b in range(a + 1, K):
            r = np.zeros(K); r[b], r[a] = 1, -1
            rows.append(r); rhs.append(offsets[a, b])
    rows.append(np.ones(K)); rhs.append(np.zeros(2))
    A, y = np.array(rows), np.array(rhs)
    t, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = A[:-1] @ t - y[:-1]
    return t, float(np.sqrt((resid ** 2).sum(1).mean()))


def excess_neighbours(g, rho, edges=EDGES, rmax=30.0, base=BASELINE):
    """Mean number of excess molecules within rmax of a molecule (model-free):
    rho * integral_0^rmax (g(r) - b) 2 pi r dr, with b the mean g over `base`."""
    g = np.atleast_2d(g)
    c = (edges[:-1] + edges[1:]) / 2
    ann = np.pi * np.diff(edges ** 2)
    sel = (c >= base[0]) & (c < base[1])
    b = (g[:, sel] * ann[sel]).sum(1) / ann[sel].sum()
    inner = edges[1:] <= rmax
    out = rho * ((g[:, inner] - b[:, None]) * ann[inner]).sum(1)
    return out if len(out) > 1 else float(out[0])


def plateau_density(xy, annulus=PLATEAU):
    """Pair-weighted point density (per nm^2) on the plateau annulus: the density
    scale implied by plateau normalization."""
    t = cKDTree(xy)
    cnt = t.count_neighbors(t, np.asarray(annulus, float))
    return float(np.diff(cnt)[0] / len(xy) / (np.pi * (annulus[1] ** 2 - annulus[0] ** 2)))


def excess_with_ci(res, rho, rng, B=400):
    """Excess neighbours within 30 nm and a block-bootstrap 95% interval."""
    nd, pr = res['null_density'], res['plateau_ratio']
    fn = lambda N, n: excess_neighbours(g_of_r(cross_g(N, n)[1], nd, pr), rho)
    val = fn(res['N'], res['n'])
    boot = bootstrap(fn, res['N'], res['n'], rng, B)
    return float(val), np.percentile(boot, [2.5, 97.5]).tolist()


def triplet_m2(xy, label, res, rng, B=400, radius=TRIPLET_RADIUS, annulus=TRIPLET_ANNULUS):
    """Connected-triplet m2 (unordered labelled partner pairs within `radius` per
    molecule) with a block-bootstrap 95% interval. Requires >= 3 labels in use."""
    k = int(label.max()) + 1
    h = baseline_excess(res['pair_g'][:k, :k])
    L = lens_integral(EDGES, h, radius)
    np.fill_diagonal(L, 0.0)
    q = connected_triplets(xy, label, k, radius, annulus, L)
    p3 = label_moments(np.bincount(label, minlength=k) / len(label))['distinct3']
    blk = res['block']
    groups = [np.flatnonzero(blk == j) for j in range(blk.max() + 1)]
    sums = np.array([q[g].sum() for g in groups]); cnts = np.array([len(g) for g in groups])
    idx = rng.integers(0, len(groups), (B, len(groups)))
    boot = sums[idx].sum(1) / cnts[idx].sum(1) / (2 * p3)
    return float(q.mean() / (2 * p3)), np.percentile(boot, [2.5, 97.5]).tolist(), float(p3)


def triplet_angles(xy, label, res, rng, B=400, radius=TRIPLET_RADIUS, annulus=TRIPLET_ANNULUS,
                   angle_edges=ANGLE_EDGES):
    """Mean per-molecule angle profile of connected three-label triplets (ordered
    neighbour pairs) with a block-bootstrap distribution."""
    k = int(label.max()) + 1
    h = baseline_excess(res['pair_g'][:k, :k])
    La = lens_angle_profile(EDGES, h, radius, angle_edges, rng)
    for b in range(k):
        La[b, b] = 0.0
    q = connected_triplet_angles(xy, label, k, radius, annulus, La, angle_edges)
    blk = res['block']
    groups = [np.flatnonzero(blk == j) for j in range(blk.max() + 1)]
    sums = np.array([q[g].sum(0) for g in groups]); cnts = np.array([len(g) for g in groups])
    idx = rng.integers(0, len(groups), (B, len(groups)))
    boot = sums[idx].sum(1) / cnts[idx].sum(1)[:, None]
    return q.mean(0), boot
