"""Cross-label statistics of the public CD20 RESI acquisitions (no clustering).

All statistics use raw round localizations. Labels are acquisition rounds
(docking sequences R1-R4). Pair counts are accumulated per spatial block so that
uncertainty can be estimated by resampling blocks.
"""
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

VIEWS = ['fig4untreated', 'edf10rep1', 'edf10rep2', 'rtx_rep1', 'rtx_rep2', 'rtx_rep3']
TREATED = {'fig4untreated': False, 'edf10rep1': False, 'edf10rep2': False,
           'rtx_rep1': True, 'rtx_rep2': True, 'rtx_rep3': True}
K = 4
LP_MAX = 8.0                      # nm, per-localization precision gate
EDGES = np.r_[np.arange(0, 60, 1.0), np.arange(60, 101, 5.0)]
PLATEAU = np.array([300.0, 600.0])  # nm, normalization annulus
PLATEAU_FRACTION = 0.05           # source subsample for the plateau counts
BLOCK = 2000.0                    # nm, spatial resampling unit


def load(cache, view, frame_fraction=1.0, lp_max=LP_MAX):
    """Concatenated localizations of the four rounds with their round label."""
    xy, lab, frame, lp = [], [], [], []
    for r in range(1, K + 1):
        d = np.load(Path(cache) / f'{view}-locs-R{r}.npz')
        keep = d['lp'] <= lp_max
        if frame_fraction < 1:
            keep &= d['frame'] < frame_fraction * (d['frame'].max() + 1)
        xy.append(np.column_stack([d['x'][keep], d['y'][keep]]).astype(float))
        lab.append(np.full(keep.sum(), r - 1))
        frame.append(d['frame'][keep]); lp.append(d['lp'][keep])
    return dict(xy=np.concatenate(xy), label=np.concatenate(lab),
                frame=np.concatenate(frame), lp=np.concatenate(lp).astype(float))


def load_resi(cache, view):
    out = []
    for r in range(1, K + 1):
        d = np.load(Path(cache) / f'{view}-resi-R{r}.npz')
        out.append(dict(xy=np.column_stack([d['x'], d['y']]), lp=d['lp']))
    return out


def blocks(xy, size=BLOCK):
    ij = np.floor((xy - xy.min(0)) / size).astype(int)
    _, idx = np.unique(ij, axis=0, return_inverse=True)
    return idx.ravel()


def block_pair_counts(xy, label, edges=EDGES, block=None, shift=None):
    """Ordered pair counts N[k, a, b, bin] with the first point in block k.

    With `shift`, every label's points are displaced by shift[label] (null
    construction: relative displacements shift[b] - shift[a]); the self-pair
    correction is then not applied.
    """
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
    """Pooled cross-label correlation per label pair, normalized by counts."""
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
    Ns, ns = N.sum(0), n.sum(0)
    _, pooled = cross_g(N, n)
    S = np.array([Ns[a, a] - ns[a] ** 2 * pooled for a in range(K)])
    return S, pooled


NULL_SMOOTH = (20.0, 100.0)       # nm, annulus defining the null pair density


def cross_annuli(xy, label, seed, shift=None, fraction=PLATEAU_FRACTION,
                 annuli=(NULL_SMOOTH, tuple(PLATEAU))):
    """Pooled cross-label pair counts per annulus from a fixed source subsample,
    normalized by sum_{a != b} n_a(subsample) n_b like `cross_g`."""
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


def null_shifts(rng, n_shift=6, lo=1500.0, hi=3000.0):
    """Labels displaced to evenly spaced directions at a common random radius, so
    every pair of labels is separated by at least 2 lo sin(pi / K)."""
    out = []
    for _ in range(n_shift):
        ang = rng.uniform(0, 2 * np.pi) + 2 * np.pi * np.arange(K) / K
        out.append(rng.uniform(lo, hi) * np.column_stack([np.cos(ang), np.sin(ang)]))
    return out


def null_normalization(xy, label, rng, n_shift=6):
    """Null pair density and plateau ratio from shifted configurations.

    Returns (null_density, plateau_ratio): the shifted-null cross pair density per
    unit annulus area over NULL_SMOOTH (pairs / sum n_a n_b / nm^2), and the
    observed/null ratio of pair counts on the plateau annulus. Shifting removes
    all cross-label correlation; the plateau ratio removes the shifted
    configuration's overlap loss.
    """
    seed = int(rng.integers(1 << 31))
    obs = cross_annuli(xy, label, seed)
    null = np.mean([cross_annuli(xy, label, seed, shift=s) for s in null_shifts(rng, n_shift)], 0)
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


def bootstrap(fn, N, n, rng, B=200):
    """Resample spatial blocks; fn(N, n) -> array."""
    nb = N.shape[0]
    out = []
    for _ in range(B):
        k = rng.integers(0, nb, nb)
        out.append(fn(N[k], n[k]))
    return np.array(out)


def registration_offsets(xy, label, radius=25.0, fraction=0.05, seed=0):
    """Mean cross-label displacement (b minus a) among close pairs, per label pair,
    from a source subsample.

    Between-molecule displacements are symmetric in expectation, so a nonzero mean
    measures a relative round registration offset.
    """
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


def excess_neighbours(g, rho, edges=EDGES, rmax=30.0, base=(60.0, 100.0)):
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
