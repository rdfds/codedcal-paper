"""Conditional covariance propagation for the SIMULATION-ONLY retrofit.

Full relative rank is a prerequisite, not a useful-precision guarantee.
The caller must supply an independently justified contrast covariance and
scientifically chosen tolerance. Missing error components are not inferred.
"""
import numpy as np
from .retrofit import t_map


def assess(c, W, kappa, contrast_covariance, *, max_sd, contrasts=None):
    """Solve centered budget-plus-blur contrasts and propagate supplied noise.

    `contrast_covariance` describes noise in c, including correlations.
    `contrasts` has one row per requested functional of all original plus
    randomized budgets; every row must sum to zero. Default: centered budgets.
    `max_sd` is an explicit SD tolerance in budget units (nm² for variance
    budgets); covariance then has units nm⁴. Passing is conditional on the
    supplied covariance, not a confidence/accuracy or physical-validity gate.
    """
    R = t_map(W, kappa)
    n = R.shape[1]
    if np.linalg.matrix_rank(R, tol=1e-10) != R.shape[0]:
        raise ValueError('not all relative budget-plus-blur components are identifiable')
    c = np.asarray(c, float)
    V = np.asarray(contrast_covariance, float)
    if c.shape != (R.shape[0],) or not np.isfinite(c).all():
        raise ValueError('one finite contrast per retained bit required')
    if V.shape != (len(c), len(c)) or not np.isfinite(V).all():
        raise ValueError('a finite square covariance for all measured contrasts is required')
    scale = max(1., float(np.max(np.abs(V))))
    if not np.allclose(V, V.T, rtol=0, atol=1e-12*scale):
        raise ValueError('contrast covariance must be symmetric')
    V = (V+V.T)/2
    if np.linalg.eigvalsh(V).min() < -1e-12*scale:
        raise ValueError('contrast covariance must be positive semidefinite')
    if not np.isscalar(max_sd) or not np.isfinite(max_sd) or max_sd <= 0:
        raise ValueError('max_sd must be an explicit finite positive tolerance')
    C = np.eye(n)-np.ones((n,n))/n if contrasts is None else np.asarray(contrasts,float)
    if (C.ndim != 2 or C.shape[1] != n or C.shape[0] == 0
            or not np.isfinite(C).all() or not np.allclose(C.sum(1),0,rtol=0,atol=1e-12)):
        raise ValueError('requested budget contrasts must be finite nonempty zero-sum rows')
    inverse = np.linalg.pinv(R)
    budget = inverse@c
    budget_covariance = inverse@V@inverse.T
    requested_covariance = C@budget_covariance@C.T
    sd = np.sqrt(np.maximum(np.diag(requested_covariance),0))
    return dict(budgets=budget, budget_covariance=budget_covariance,
                contrast_estimates=C@budget, contrast_covariance=requested_covariance,
                contrast_sd=sd, precision_rule_passed=bool(np.all(sd <= max_sd)),
                max_sd=float(max_sd), estimand='relative budget PLUS per-plane blur',
                uncertainty_scope='conditional propagation of supplied covariance; '
                    'omitted noise, model bias and physical validity are not certified',
                status='SIMULATION ONLY; no physical randomized MERFISH readout validated')


def abundance_path(W, start_counts, end_counts):
    """Find certified rank loss on one straight abundance path.

    A change in inertia of H=W.T diag(k) W-diag(W.T k) proves a zero
    eigenvalue between endpoints. Since P1=1 and R=[-P,1], rank(R)=rank(P)
    and P=D^-1 H/(h-1). All bit abundances must be positive throughout.
    This diagnoses a failure missed by checking only observed profiles.
    Finding no crossing DOES NOT certify path or convex-hull robustness:
    eigenvalues can cross twice, and other paths remain unexamined.
    """
    from scipy.optimize import brentq
    W = np.asarray(W, float)
    if (W.ndim != 2 or not np.isfinite(W).all() or not np.isin(W,[0,1]).all()
            or W.shape[0] == 0 or W.shape[1] < 2):
        raise ValueError('a nonempty binary codebook is required')
    weights = np.unique(W.sum(1))
    if len(weights) != 1 or weights[0] < 2:
        raise ValueError('one constant code weight of at least two is required')
    start = np.asarray(start_counts,float); end = np.asarray(end_counts,float)
    for counts in (start,end):
        if (counts.shape != (len(W),) or not np.isfinite(counts).all()
                or (counts < 0).any() or not np.isfinite(counts@W).all()
                or ((counts@W) <= 0).any()):
            raise ValueError('finite nonnegative counts with positive abundance in every bit required')
    def eigenvalues(fraction):
        counts=(1-fraction)*start+fraction*end
        total=counts@W
        H=W.T@(counts[:,None]*W)-np.diag(total)
        return np.linalg.eigvalsh(H)
    left=eigenvalues(0.);right=eigenvalues(1.);roots=[]
    for fraction, values in [(0.,left),(1.,right)]:
        singular=np.linalg.svd(t_map(W,(1-fraction)*start+fraction*end),compute_uv=False)
        if singular[-1] < 1e-10:
            roots.append(dict(eigen_index=int(np.argmin(np.abs(values))),fraction=fraction,
                              smallest_singular=float(singular[-1])))
    changes=(np.sign(left)!=np.sign(right)) & (left!=0) & (right!=0)
    for index in np.flatnonzero(changes):
        fraction=brentq(lambda t:eigenvalues(t)[index],0.,1.,xtol=1e-14)
        counts=(1-fraction)*start+fraction*end
        singular=np.linalg.svd(t_map(W,counts),compute_uv=False)
        if not any(abs(fraction-r['fraction'])<1e-10 for r in roots):
            roots.append(dict(eigen_index=int(index),fraction=float(fraction),
                              smallest_singular=float(singular[-1])))
    return dict(endpoint_inertias=[int(sum(left<0)),int(sum(right<0))],
                rank_loss_found=bool(roots),crossings=roots,
                robustness_certified=False,
                scope='sufficient diagnostic of rank loss; absence of crossings is not certification',
                status='SIMULATION ONLY; no physical randomized MERFISH readout validated')
