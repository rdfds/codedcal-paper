"""Nuclear pore complexes as a stoichiometric ground truth for label statistics.

Each Nup96 NPC holds n = 32 copies. A copy is observed with label a with
probability q_a = p pi_a, independently across copies (labelling and detection
efficiency p, observed label frequencies pi). For per-pore label counts N_a:

    E[sum_{a<b} N_a N_b]       = n(n-1)      p^2 e2(pi)
    E[sum_{a<b<c} N_a N_b N_c] = n(n-1)(n-2) p^3 e3(pi)

(e_k: elementary symmetric polynomials). Cross-label pairs and three-distinct-
label triples never pair a molecule with itself, so they need no resolution of
same-label neighbours. Each gives an estimate of p for known n; together they
also give n and p jointly.

Corner level: within radius R of a molecule in the same pore, the cross-label
count and the three-distinct-label pair count estimate p a1 and p^2 a2, where a1
and a2 are the expected numbers of other copies and of copy pairs within R
(from the pore geometry), so m2 / n_ex = p a2 / a1.
"""
from math import comb
import numpy as np
from scipy.optimize import brentq
from scipy.spatial import cKDTree

N_COPIES = 32
RING_RADIUS = 51.0        # nm; cryo-EM neighbour spacings 39 nm (adjacent corners), ~101 nm across
PAIR_SPACING = 12.0       # nm, two copies per ring per corner
RING_OFFSET = 3.0         # nm, radial offset between cytoplasmic and nuclear copies
K = 4


def elementary(pi, k):
    e = np.zeros(k + 1); e[0] = 1.0
    for p in pi:
        e[1:] = e[1:] + p * e[:-1].copy()
    return e[k]


def corner_template():
    """Top-view positions of the four Nup96 copies of one corner, relative to the
    corner centre (tangential = x, radial = y)."""
    h = PAIR_SPACING / 2
    return np.array([[-h, RING_OFFSET], [h, RING_OFFSET], [-h, -RING_OFFSET], [h, -RING_OFFSET]])


def pore_positions(rng, rotation=None):
    """32 copy positions of one pore centred at the origin."""
    rot = rng.uniform(0, 2 * np.pi) if rotation is None else rotation
    pts = []
    t = corner_template()
    for c in range(8):
        a = rot + 2 * np.pi * c / 8
        tang, rad = np.array([-np.sin(a), np.cos(a)]), np.array([np.cos(a), np.sin(a)])
        centre = RING_RADIUS * rad
        pts.append(centre + t[:, :1] * tang + t[:, 1:] * rad)
    return np.concatenate(pts)


def simulate(rng, n_pores, p, pi, link_sd=4.0, loc_sd=0.7, merge_nm=0.0, spacing=400.0):
    """Observed molecules of n_pores pores on a grid: positions, labels, pore ids.

    Copies are observed with probability p, labelled by pi; same-label molecules
    of one pore closer than merge_nm are merged (cluster-merging artefact)."""
    side = int(np.ceil(np.sqrt(n_pores)))
    xy, lab, pid = [], [], []
    for k in range(n_pores):
        centre = spacing * np.array([k % side, k // side], float)
        pos = pore_positions(rng) + centre
        keep = rng.random(N_COPIES) < p
        pos = pos[keep] + rng.normal(0, link_sd, (keep.sum(), 2)) + rng.normal(0, loc_sd, (keep.sum(), 2))
        lb = rng.choice(K, keep.sum(), p=pi)
        if merge_nm > 0:
            pos, lb = _merge(pos, lb, merge_nm)
        xy.append(pos); lab.append(lb); pid.append(np.full(len(pos), k))
    return np.concatenate(xy), np.concatenate(lab), np.concatenate(pid)


def _merge(pos, lb, d):
    out_p, out_l = [], []
    for a in np.unique(lb):
        P = pos[lb == a]
        alive = np.ones(len(P), bool)
        for i in range(len(P)):
            if not alive[i]:
                continue
            close = np.flatnonzero(alive & (np.hypot(*(P - P[i]).T) < d))
            out_p.append(P[close].mean(0)); out_l.append(a)
            alive[close] = False
    return np.array(out_p).reshape(-1, 2), np.array(out_l, int)


def geometry_factors(rng, radius, link_sd=4.0, loc_sd=0.7, n_pores=400):
    """a1, a2: expected other copies and copy pairs within `radius` of a copy in
    the same pore, with every copy observed."""
    a1, a2 = [], []
    for _ in range(n_pores):
        pos = pore_positions(rng) + rng.normal(0, link_sd, (N_COPIES, 2)) + rng.normal(0, loc_sd, (N_COPIES, 2))
        c = cKDTree(pos).query_ball_point(pos, radius, return_length=True) - 1
        a1.append(c.mean()); a2.append((c * (c - 1) / 2).mean())
    return float(np.mean(a1)), float(np.mean(a2))


def pore_counts(lab, pid):
    """Per-pore label counts N[k, a] (pores with at least one molecule)."""
    u, inv = np.unique(pid, return_inverse=True)
    N = np.zeros((len(u), K))
    np.add.at(N, (inv, lab), 1)
    return N


def factorial_moments(N, pi):
    """Mean cross-label pairs / e2 and three-distinct-label triples / e3 per pore."""
    X2 = sum(N[:, a] * N[:, b] for a in range(K) for b in range(a + 1, K))
    X3 = sum(N[:, a] * N[:, b] * N[:, c] for a in range(K) for b in range(a + 1, K) for c in range(b + 1, K))
    return X2.mean() / elementary(pi, 2), X3.mean() / elementary(pi, 3)


def stoichiometry(N, pi, n=N_COPIES):
    """p from counts, pairs and triples for known n, and (n, p) jointly from pairs+triples."""
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


def corner_statistics(xy, lab, pid, pi, radius):
    """n_ex (labelled partners within radius, from cross-label counts) and m2
    (labelled partner pairs, from three-distinct-label pairs), same pore only.
    Returns per-molecule contributions for resampling by pore."""
    pi = np.asarray(pi)
    s2, s3 = (pi ** 2).sum(), (pi ** 3).sum()
    p3 = 1 - 3 * s2 + 2 * s3
    tree = cKDTree(xy)
    nb = tree.query_ball_point(xy, radius)
    c = np.zeros(len(xy)); t = np.zeros(len(xy))
    for i, js in enumerate(nb):
        js = np.asarray(js, int)
        js = js[(pid[js] == pid[i]) & (lab[js] != lab[i])]
        if not len(js):
            continue
        cnt = np.bincount(lab[js], minlength=K).astype(float)
        c[i] = cnt.sum()
        t[i] = cnt.sum() ** 2 - (cnt ** 2).sum()
    return c / (1 - s2), t / (2 * p3)


def resample_pores(values_per_mol, pid, rng, B=400):
    """Bootstrap over pores of per-molecule means; values_per_mol: list of arrays."""
    u, inv = np.unique(pid, return_inverse=True)
    sums = [np.bincount(inv, weights=v, minlength=len(u)) for v in values_per_mol]
    cnt = np.bincount(inv, minlength=len(u)).astype(float)
    out = []
    for _ in range(B):
        k = rng.integers(0, len(u), len(u))
        out.append([s[k].sum() / cnt[k].sum() for s in sums])
    return np.array(out)
