"""SIMULATION ONLY: abundance-weighted randomized-readout contrasts.

Constant-weight codebooks permit exact mixture cancellation under stable
weights and the additive second-moment model. Recovering all relative
budget-PLUS-blur components additionally requires full row rank of t_map.
A randomized round does not guarantee that rank for every codebook.

Separating displacement budgets from arbitrary plane blur requires an
explicit observation model. assess_blur_model supplies numerical witnesses
for such models; a single replicate is not sufficient in general. No
empirically randomized MERFISH acquisition is validated by this module.
"""
import itertools

import numpy as np


def edge_structure(W):
    """Bit pairs with nonzero gene overlap; A = (edges x genes) incidence."""
    nb = W.shape[1]
    edges_all = list(itertools.combinations(range(nb), 2))
    A_all = np.array([(W[:, i] * W[:, j]) for i, j in edges_all])
    keep = A_all.sum(1) > 0
    return [e for e, k in zip(edges_all, keep) if k], A_all[keep]


def contrast_weights(W, kappa):
    """w_ij = K_ij / ((h-1) K_i): abundance-weighted contrast weights.

    `kappa` is per-gene abundance (decoded counts are the intended input).
    The 3 is C(weight-1, 1) for a weight-4 book: each bit i sits in 3 pairs
    per codeword, so the weights sum to one over the pairs containing i.
    """
    W = np.asarray(W)
    kappa = np.asarray(kappa, float)
    if (W.ndim != 2 or min(W.shape) == 0 or not np.isin(W, [0, 1]).all()
            or kappa.shape != (len(W),) or not np.isfinite(kappa).all()
            or np.any(kappa < 0)):
        raise ValueError('binary codebook and finite nonnegative target abundances required')
    W = W.astype(int)
    if np.any(W.sum(1) < 2):
        raise ValueError('each codeword must have weight at least two')
    if np.any(W.T @ kappa <= 0):
        raise ValueError('every retained bit must have positive abundance')
    edges, A = edge_structure(W)
    nb = W.shape[1]
    K_i = W.T.astype(float) @ kappa
    K_ij = A @ kappa
    weights = np.unique(W.sum(1))
    weights = weights[weights > 0]
    if len(weights) != 1:
        raise ValueError(
            'the contrast weights assume a CONSTANT-WEIGHT codebook (each bit '
            f'sits in the same number of pairs per codeword); this book has '
            f'weights {weights.tolist()}. Drop non-coding rows or supply a '
            'constant-weight subset.')
    per_bit = int(weights[0]) - 1         # pairs per codeword containing bit i
    Wt = np.zeros((nb, len(edges)))
    for e, (i, j) in enumerate(edges):
        if K_i[i] > 0:
            Wt[i, e] = K_ij[e] / (per_bit * K_i[i])
        if K_i[j] > 0:
            Wt[j, e] = K_ij[e] / (per_bit * K_i[j])
    return Wt, edges


def budget_contrasts(m_pairs, m_random, W, kappa):
    """Per-round budget contrasts from edge statistics.

    m_pairs   : edge statistics m[i,j], ordered as `edge_structure(W)` edges.
    m_random  : per-bit statistics m[i,R] against the randomized round R.
    Returns c (n_bits,), the mixture-free contrasts.
    """
    Wt, edges = contrast_weights(W, kappa)
    m_pairs = np.asarray(m_pairs, float)
    m_random = np.asarray(m_random, float)
    if m_pairs.shape[0] != len(edges):
        raise ValueError(f'm_pairs has {m_pairs.shape[0]} entries, '
                         f'expected {len(edges)} edges')
    return m_random - Wt @ m_pairs


def t_map(W, kappa):
    """The linear map R with c = R @ [T_0..T_{nb-1}, T_R].

    Because sum_j w_ij = 1, the contrast reduces to
        c_i = T_R - sum_j w_ij T_j,
    i.e. c is NOT T_i directly but the randomized round's budget minus an
    abundance-weighted average of the partner planes'. Budget contrasts are
    recovered by inverting this map (`solve_budgets`), which is why the map's
    rank is the identifiability statement.
    """
    Wt, edges = contrast_weights(W, kappa)
    nb = np.asarray(W).shape[1]
    R = np.zeros((nb, nb + 1))
    R[:, nb] = 1.
    for e, (i, j) in enumerate(edges):
        R[i, j] -= Wt[i, e]
        R[j, i] -= Wt[j, e]
    return R


def solve_budgets(c, W, kappa):
    """Recover per-plane budget CONTRASTS from the measured contrasts c.

    Returns T_hat (n_bits + 1,) centred so its mean is zero: the global level
    is not identified. Full contrast rank is checked before inversion;
    otherwise the solver refuses to invent unidentifiable differences. The last entry is the randomized round itself.
    """
    R = t_map(W, kappa)
    if np.linalg.matrix_rank(R, tol=1e-10) != R.shape[0]:
        raise ValueError('not all budget-plus-blur contrasts are identifiable for this book')
    c = np.asarray(c, float)
    if c.shape != (R.shape[0],) or not np.isfinite(c).all():
        raise ValueError('one finite contrast per retained bit required')
    T, *_ = np.linalg.lstsq(R, c, rcond=None)
    return T - T.mean()


def identified_map(W, kappa):
    """The map T -> c, its rank, and its null space.

    rank == n_bits with a constant null vector means every per-round budget
    contrast is identified and only the global level is unrecoverable.
    """
    R = t_map(W, kappa)
    nb = np.asarray(W).shape[1]
    rank = int(np.linalg.matrix_rank(R, tol=1e-10))
    null = np.linalg.svd(R)[2][-1]
    return dict(rank=rank, n_budgets=nb + 1,
                null_vector=null,
                null_is_global_constant=float(np.abs(null - null.mean()).max()),
                all_contrasts_identified=bool(rank == nb))


def check_mixture_cancellation(W, kappa, sigma2):
    """Verify the contrast kills an arbitrary per-gene sigma^2 mixture.

    Returns the max absolute residual; it should be at machine precision for
    ANY sigma2, which is the whole point of the construction.
    """
    W = np.asarray(W, int)
    kappa = np.asarray(kappa, float)
    sigma2 = np.asarray(sigma2, float)
    edges, A = edge_structure(W)
    Wt, _ = contrast_weights(W, kappa)
    K_i = W.T.astype(float) @ kappa
    K_ij = A @ kappa
    m_d2 = (A @ (kappa * sigma2)) / np.maximum(K_ij, 1e-300)
    m_d1 = (W.T.astype(float) @ (kappa * sigma2)) / np.maximum(K_i, 1e-300)
    return float(np.abs(m_d1 - Wt @ m_d2).max())


def contrast_noise(W, n_spots_per_gene, kappa=None, sigma2=None, n_trials=200,
                   seed=0):
    """Finite-N noise on the contrast, from binomial randomization alone.

    The extra readout hits a Binomial(N, 1/2) share of each gene's spots, so
    the realized abundance behind m[i,R] fluctuates. Returns the sd of the
    contrast in the same units as sigma2 (nm^2 with the defaults), which
    falls as 1/sqrt(N) — a noise, not a floor.
    """
    rng = np.random.default_rng(seed)
    W = np.asarray(W, int)
    G = W.shape[0]
    edges, A = edge_structure(W)
    out = []
    for _ in range(n_trials):
        kap = rng.lognormal(0, 1.3, G) if kappa is None else np.asarray(kappa, float)
        s2 = rng.uniform(400., 3025., G) if sigma2 is None else np.asarray(sigma2, float)
        kap_r = kap * (rng.binomial(n_spots_per_gene, 0.5, G)
                       / (n_spots_per_gene / 2.))
        K_ij = A @ kap
        m_d2 = (A @ (kap * s2)) / np.maximum(K_ij, 1e-300)
        K_i_r = W.T.astype(float) @ kap_r
        m_d1 = (W.T.astype(float) @ (kap_r * s2)) / np.maximum(K_i_r, 1e-300)
        out.append(budget_contrasts(m_d2, m_d1, W, kap))
    out = np.array(out)
    return dict(n_spots_per_gene=int(n_spots_per_gene),
                sd=float(np.std(out, axis=0).mean()),
                max_abs_mean=float(np.abs(out.mean(0)).max()))


def report(W, kappa=None, n_spots=(100, 1000), seed=0):
    """Design-time summary: what one randomized round would buy this book."""
    rng = np.random.default_rng(seed)
    G = np.asarray(W).shape[0]
    kap = rng.lognormal(0, 1.3, G) if kappa is None else np.asarray(kappa, float)
    s2 = rng.uniform(400., 3025., G)
    return dict(
        n_genes=int(G), n_bits=int(np.asarray(W).shape[1]),
        mixture_cancellation_max_abs=check_mixture_cancellation(W, kap, s2),
        identified=identified_map(W, kap),
        noise={f'N_{n}': contrast_noise(W, n, seed=seed + n) for n in n_spots},
        estimand='budget PLUS per-plane blur; separation requires a specified '
                 'augmented model, not merely one replicated plane',
        blur_separation=assess_blur_model(W, kap),
        status='SIMULATION ONLY - this project has no experimentally validated '
               'per-molecule randomized MERFISH readout')


def assess_blur_model(W, kappa, blur_basis=None, blur_measurements=None):
    """Assess relative T in c=R(T+B alpha), v=H B alpha.

    B defaults to independent blur for every original and randomized plane.
    H declares *independent blur observations*, NOT repeat-image labels.
    Their physical validity must be justified separately. B=ones assumes
    common blur; it does not infer common blur from the data.
    """
    from .identifiability import assess, relative_contrasts
    R = t_map(W, kappa)
    n = R.shape[1]
    B = np.eye(n) if blur_basis is None else np.asarray(blur_basis, float)
    H = np.empty((0, n)) if blur_measurements is None else np.asarray(blur_measurements, float)
    if (B.ndim != 2 or B.shape[0] != n or H.ndim != 2 or H.shape[1] != n
            or not np.isfinite(B).all() or not np.isfinite(H).all()):
        raise ValueError('blur basis and observation dimensions must match all budgets')
    A = np.block([[R, R @ B], [np.zeros((len(H), n)), H @ B]])
    C = np.column_stack([relative_contrasts(n), np.zeros((n-1, B.shape[1]))])
    out = assess(A, C)
    out.update(model='c=R(T+B alpha); v=H B alpha',
               parameter_order='all T, then blur-basis coefficients alpha',
               physical_blur_observations_verified=False,
               blur_basis=B.tolist(), blur_measurements=H.tolist())
    return out
