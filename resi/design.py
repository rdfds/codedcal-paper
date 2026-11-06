"""Design of stochastically multiplexed acquisitions: number and length of rounds.

An m-molecule configuration contributes to an m-th order cross-label statistic
only when all m molecules are labelled, detected and carry distinct labels:

    P(contributes) = (efficiency * p_det)^m * D_m(pi),
    D_m(pi) = m! e_m(pi)      (e_m: elementary symmetric polynomial),

which for K equiprobable labels is prod_{j<m} (1 - j/K). Splitting a fixed total
acquisition into K rounds shortens each round to 1/K of the total, which lowers
the per-molecule detection probability p_det through the minimum localization
count required to call a molecule. The information available for an m-th order
question therefore scales as (p_det(T/K))^m D_m(K).
"""
from math import factorial
import numpy as np


def distinct_probability(m, pi):
    """Probability that m molecules with independent labels ~ pi carry distinct labels."""
    pi = np.asarray(pi, float) / np.sum(pi)
    e = np.zeros(m + 1)
    e[0] = 1.0
    for p in pi:
        e[1:] = e[1:] + p * e[:-1].copy()
    return factorial(m) * e[m]


def uniform_distinct(m, K):
    return float(np.prod([1 - j / K for j in range(m)])) if K >= m else 0.0


def relative_information(m, K, phi, p_det, K_ref=4, phi_ref=1.0):
    """m-th order information of K rounds of relative length phi, relative to the
    reference design; p_det(phi) is the detection probability vs round length."""
    return (p_det(phi) / p_det(phi_ref)) ** m * uniform_distinct(m, K) / uniform_distinct(m, K_ref)


def fixed_budget(m, K_values, p_det, K_ref=4):
    """Relative information when a fixed total time (that of K_ref full rounds) is
    split into K equal rounds: each round has relative length phi = K_ref / K."""
    return {int(K): relative_information(m, K, K_ref / K, p_det, K_ref) for K in K_values}


def interpolated(phis, pdets):
    """p_det(phi) by log-linear interpolation of measured points, flat beyond them."""
    lp, pd = np.log(np.asarray(phis, float)), np.asarray(pdets, float)
    order = np.argsort(lp)
    lp, pd = lp[order], pd[order]
    return lambda phi: float(np.interp(np.log(phi), lp, pd))
