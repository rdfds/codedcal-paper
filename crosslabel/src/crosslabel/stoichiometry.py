"""Copy number and labelling efficiency of discrete assemblies from label statistics.

Each assembly (group, e.g. one nuclear pore) holds n copies. A copy is observed
with label a with probability q_a = p pi_a, independently across copies
(labelling-and-detection efficiency p, observed label frequencies pi). For
per-group label counts N_a:

    E[sum_{a<b} N_a N_b]       = n(n-1)      p^2 e2(pi)
    E[sum_{a<b<c} N_a N_b N_c] = n(n-1)(n-2) p^3 e3(pi)

(e_k: elementary symmetric polynomials). Cross-label pairs and three-distinct-
label triples never pair a molecule with itself, so they need no resolution of
same-label neighbours. Each gives an estimate of p for known n; together they
give n and p jointly.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.spatial import cKDTree


def elementary(pi, k):
    """Elementary symmetric polynomial e_k of the entries of pi."""
    e = np.zeros(k + 1); e[0] = 1.0
    for p in pi:
        e[1:] = e[1:] + p * e[:-1].copy()
    return e[k]


def group_counts(label, group, n_labels=None):
    """Per-group label counts N[g, a] (groups with at least one point)."""
    K = int(label.max()) + 1 if n_labels is None else int(n_labels)
    u, inv = np.unique(group, return_inverse=True)
    N = np.zeros((len(u), K))
    np.add.at(N, (inv, label), 1)
    return N


def factorial_moments(N, pi):
    """Mean cross-label pairs / e2 and three-distinct-label triples / e3 per group:
    estimates of n(n-1) p^2 and n(n-1)(n-2) p^3."""
    K = N.shape[1]
    X2 = sum(N[:, a] * N[:, b] for a in range(K) for b in range(a + 1, K))
    X3 = sum(N[:, a] * N[:, b] * N[:, c] for a in range(K) for b in range(a + 1, K) for c in range(b + 1, K))
    return X2.mean() / elementary(pi, 2), X3.mean() / elementary(pi, 3)


def stoichiometry(N, pi, n):
    """Efficiency from counts, cross-label pairs and three-label triples for known
    copy number n, and (n, p) jointly from pairs and triples."""
    f2, f3 = factorial_moments(N, pi)
    p1 = N.sum(1).mean() / n
    p2 = np.sqrt(f2 / (n * (n - 1)))
    p3 = (f3 / (n * (n - 1) * (n - 2))) ** (1 / 3)
    r = f3 / f2                       # = (n - 2) p
    try:
        n_joint = brentq(lambda m: m * (m - 1) * (r / (m - 2)) ** 2 - f2, 3.0 + 1e-6, 1e4)
        p_joint = r / (n_joint - 2)
    except ValueError:
        n_joint, p_joint = np.nan, np.nan
    return dict(p_counts=p1, p_pairs=p2, p_triples=p3, n_joint=n_joint, p_joint=p_joint)


def group_neighbour_statistics(xy, label, group, pi, radius, n_labels=None):
    """Within-group neighbours of other labels within `radius`, per molecule:
    returns (n_ex, m2) contributions, the labelled partners and labelled partner
    pairs implied by cross-label and three-distinct-label counts."""
    K = int(label.max()) + 1 if n_labels is None else int(n_labels)
    pi = np.asarray(pi)
    s2, s3 = (pi ** 2).sum(), (pi ** 3).sum()
    p3 = 1 - 3 * s2 + 2 * s3
    tree = cKDTree(xy)
    nb = tree.query_ball_point(xy, radius)
    c = np.zeros(len(xy)); t = np.zeros(len(xy))
    for i, js in enumerate(nb):
        js = np.asarray(js, int)
        js = js[(group[js] == group[i]) & (label[js] != label[i])]
        if not len(js):
            continue
        cnt = np.bincount(label[js], minlength=K).astype(float)
        c[i] = cnt.sum()
        t[i] = cnt.sum() ** 2 - (cnt ** 2).sum()
    return c / (1 - s2), t / (2 * p3)
