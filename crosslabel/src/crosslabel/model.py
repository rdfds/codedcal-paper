"""Cross-label pair correlation of a monomer/dimer mixture under complete spatial randomness.

Two labelled subunits of a dimer are separated by d (random orientation). Each
localization adds the label linkage offset (per-axis sd sigma_link, fixed per
molecule) and its own localization error. The cross-label displacement between
localizations of the two partners is then d*u + N(0, s^2 I) with
s^2 = 2 sigma_link^2 + lp_a^2 + lp_b^2, a 2D Rice law with per-area density

    h(r) = exp(-(r^2 + d^2) / (2 s^2)) I0(r d / s^2) / (2 pi s^2).

Different labels imply different molecules, so the cross-label correlation is

    g(r) = 1 + (f_d / rho) h(r),

where f_d is the fraction of subunits in dimers and rho the subunit density
(labelled density / labelling efficiency). Neither the labelling efficiency nor
per-molecule localization counts enter g beyond rho, provided counts are
independent of organization.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.special import i0e


def rice_area_density(r, d, s):
    z = r * d / s ** 2
    # i0e(z) = exp(-z) I0(z), so exp(-(r^2+d^2)/2s^2) I0(z) = exp(-(r-d)^2/2s^2) i0e(z)
    return np.exp(-(r - d) ** 2 / (2 * s ** 2)) * i0e(z) / (2 * np.pi * s ** 2)


def bin_average(fun, edges, n=8):
    """Area-weighted average of fun(r) over each annular bin."""
    lo, hi = edges[:-1], edges[1:]
    t = (np.arange(n) + 0.5) / n
    r = lo[:, None] + (hi - lo)[:, None] * t
    w = r
    return (fun(r) * w).sum(1) / w.sum(1)


def dimer_g(edges, amplitude, d, sigma_link, lp_sq_sum, baseline=1.0):
    s = np.sqrt(2 * sigma_link ** 2 + lp_sq_sum)
    return baseline + amplitude * bin_average(lambda r: rice_area_density(r, d, s), edges)


def fit_dimer(edges, g, se, lp_sq_sum, rmax=100.0, free_baseline=True):
    """Weighted least squares fit of the monomer/dimer cross-label model.

    Returns amplitude (= f_d / rho, nm^2), d, sigma_link (nm) and baseline.
    """
    c = (edges[:-1] + edges[1:]) / 2
    idx = np.flatnonzero((c < rmax) & np.isfinite(g) & (se > 0))

    def model(p):
        amp, d, sl = p[:3]
        base = p[3] if free_baseline else 1.0
        full = dimer_g(edges, amp, d, sl, lp_sq_sum, base)
        return full[idx]

    def resid(p):
        return (model(p) - g[idx]) / se[idx]

    p0 = [200.0, 13.0, 4.0] + ([1.0] if free_baseline else [])
    lb = [0, 3.0, 0.5] + ([0.5] if free_baseline else [])
    ub = [1e5, 30.0, 15.0] + ([2.0] if free_baseline else [])
    fit = least_squares(resid, p0, bounds=(lb, ub))
    J = fit.jac
    try:
        cov = np.linalg.inv(J.T @ J) * max(1.0, (fit.fun ** 2).sum() / max(len(idx) - len(p0), 1))
    except np.linalg.LinAlgError:
        cov = np.full((len(p0), len(p0)), np.nan)
    names = ['amplitude_nm2', 'd_nm', 'sigma_link_nm'] + (['baseline'] if free_baseline else [])
    return dict(params=dict(zip(names, fit.x.tolist())),
                se=dict(zip(names, np.sqrt(np.diag(cov)).tolist())),
                chi2=float((fit.fun ** 2).sum()), dof=int(len(idx) - len(p0)),
                fitted=model(fit.x), r=c[idx])


def bootstrap_fit(edges, gboot, se, lp_sq_sum, **kw):
    """Refit every bootstrap replicate; returns an array of parameter vectors."""
    rows = []
    for gb in gboot:
        f = fit_dimer(edges, gb, se, lp_sq_sum, **kw)
        rows.append(list(f['params'].values()))
    return np.array(rows)
