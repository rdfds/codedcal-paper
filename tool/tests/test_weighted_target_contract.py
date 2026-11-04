"""Exact population controls for optional event weights, without sampled data."""
import itertools
from fractions import Fraction

import numpy as np

from codedcal.smlm import aggregate_events, components


def population_moments(dependent_weights=False):
    # One stationary source, iid errors in {-1,0,1}, all four permanent codes.
    # Enumeration averages the finite law exactly; these are not biological
    # replicates or confidence-coverage simulations.
    states = np.array(list(itertools.product([-1., 0., 1.], repeat=4)))
    x = np.repeat(states, 4, axis=0)
    group = np.repeat(np.arange(len(x)), 4)
    code = np.repeat(np.tile(np.arange(4), len(states)), 4)
    panel = np.tile(np.arange(4), len(x))
    if dependent_weights:
        weight = np.where(np.arange(4)[None, :] >= 2,
                          np.where(x != 0, 2., 1.), 1.)
    else:
        weight = np.tile([1., 3., 2., 5.], (len(x), 1))
    z, n, _ = aggregate_events(group, code, panel, x.ravel(), weight.ravel(),
                               q=[np.pi / 3])
    return components(z[..., 0], n).mean(axis=0)[:, :2]


def test_fixed_unequal_weights_preserve_stationary_population_target():
    # E exp(i*q*error)=2/3. For counts 1,3,2,5 the two reference
    # moments are 14/3 and 80/9, and both cross moments are identical.
    expected = np.array([[Fraction(14, 3), Fraction(14, 3)],
                         [Fraction(80, 9), Fraction(80, 9)]], dtype=float)
    np.testing.assert_allclose(population_moments(), expected, rtol=0, atol=2e-13)


def test_positive_error_dependent_weights_do_not_guarantee_invariance():
    # Zero physical motion and identical unweighted errors in every panel.
    # Only window B upweights nonzero errors. Effective weighted kernels
    # therefore change; arithmetic C/R=48/47 despite the null motion target.
    # Retain the stronger geometric alternative: its ratio is 1 in this
    # particular one-source example. That is not a general safety guarantee.
    expected = np.array([[Fraction(2, 3), Fraction(20, 27)],
                         [Fraction(10, 9), 1]], dtype=float)
    got = population_moments(dependent_weights=True)
    np.testing.assert_allclose(got, expected, rtol=0, atol=2e-13)
    arithmetic = got[:, 0].mean() / got[:, 1].mean()
    geometric = np.sqrt(np.prod(got[:, 0] / got[:, 1]))
    np.testing.assert_allclose(arithmetic, float(Fraction(48, 47)), rtol=0, atol=2e-13)
    np.testing.assert_allclose(geometric, 1., rtol=0, atol=2e-13)
