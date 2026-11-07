"""Oracle quantities (molecule identities known), used only to verify the exact
label identities on simulated data."""
import numpy as np
from scipy.spatial import cKDTree

from .labels import label_moments

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
