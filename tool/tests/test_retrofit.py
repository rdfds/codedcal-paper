"""Tests for the one-randomized-round retrofit estimator.

The numbers checked here are the ones recorded in
`outputs/code-design-sep16/d1_retrofit_results.json`; the point of the module
is that the abundance-weighted contrast cancels the per-gene sigma^2 mixture
IDENTICALLY, for any abundance and any mixture, so the tests assert machine
precision rather than a tolerance band.
"""
import numpy as np
import pytest

from codedcal import certificate, retrofit


@pytest.fixture(scope='module')
def W():
    return certificate.rm24_shell()


def test_shell_shape(W):
    assert W.shape == (140, 16)
    assert set(W.sum(1)) == {4}


def test_mixture_cancels_identically(W):
    """The contrast must kill ANY per-gene mixture, not a typical one."""
    rng = np.random.default_rng(0)
    for _ in range(25):
        kappa = rng.lognormal(0, rng.uniform(0.1, 2.5), len(W))
        sigma2 = rng.uniform(1, 1e4, len(W))
        assert retrofit.check_mixture_cancellation(W, kappa, sigma2) < 1e-9


def test_mixture_cancels_for_degenerate_abundance(W):
    """Uniform abundance and near-singular abundance are both fine."""
    n = len(W)
    assert retrofit.check_mixture_cancellation(W, np.ones(n),
                                               np.linspace(400, 3025, n)) < 1e-9
    kappa = np.full(n, 1e-6)
    kappa[0] = 1.0
    assert retrofit.check_mixture_cancellation(W, kappa,
                                               np.linspace(400, 3025, n)) < 1e-9


def test_all_contrasts_identified_null_is_constant(W):
    rng = np.random.default_rng(1)
    out = retrofit.identified_map(W, rng.lognormal(0, 1.3, len(W)))
    assert out['all_contrasts_identified'] is True
    assert out['rank'] == W.shape[1]
    assert out['n_budgets'] == W.shape[1] + 1
    # the unrecoverable direction is the GLOBAL LEVEL and nothing else
    assert out['null_is_global_constant'] < 1e-9


def test_contrast_noise_falls_as_one_over_sqrt_n(W):
    """Finite-N randomization noise is a noise, not a floor."""
    a = retrofit.contrast_noise(W, 100, n_trials=120, seed=3)
    b = retrofit.contrast_noise(W, 1000, n_trials=120, seed=3)
    assert b['sd'] < a['sd']
    ratio = a['sd'] / b['sd']
    assert 2.0 < ratio < 5.0          # sqrt(10) = 3.16 expected


def test_budget_contrasts_recovers_injected_contrasts(W):
    """With the mixture present, injected per-plane budgets are recovered up
    to the global constant."""
    rng = np.random.default_rng(7)
    G, nb = W.shape
    kappa = rng.lognormal(0, 1.3, G)
    sigma2 = rng.uniform(400., 3025., G)
    T = rng.uniform(0., 300., nb)                 # injected per-plane budgets
    T_R = 50.0                                    # randomized round's own budget
    edges, A = retrofit.edge_structure(W)
    K_ij = A @ kappa
    K_i = W.T.astype(float) @ kappa
    # edge statistics carry budget sum + the gene mixture
    m_pairs = ((A @ (kappa * sigma2)) / K_ij
               + np.array([T[i] + T[j] for i, j in edges]))
    m_random = (W.T.astype(float) @ (kappa * sigma2)) / K_i + T + T_R
    c = retrofit.budget_contrasts(m_pairs, m_random, W, kappa)
    # c is T_R minus a weighted average of partners, so the budgets are
    # recovered by INVERTING the map, not by reading c directly.
    T_hat = retrofit.solve_budgets(c, W, kappa)
    truth = np.append(T, T_R)
    truth = truth - truth.mean()
    r = np.corrcoef(T_hat, truth)[0, 1]
    assert r > 0.99, f'recovery correlation {r:.3f}'
    slope = np.polyfit(truth, T_hat, 1)[0]
    assert 0.9 < slope < 1.1, f'recovery slope {slope:.3f}'
    assert np.abs(T_hat - truth).max() < 1e-6, 'recovery should be exact'


def test_solve_budgets_is_exact_under_any_mixture(W):
    """Recovery must not depend on the per-gene mixture at all."""
    rng = np.random.default_rng(11)
    G, nb = W.shape
    kappa = rng.lognormal(0, 2.0, G)
    T = rng.uniform(0., 500., nb)
    T_R = 123.4
    edges, A = retrofit.edge_structure(W)
    K_ij, K_i = A @ kappa, W.T.astype(float) @ kappa
    for _ in range(5):
        sigma2 = rng.uniform(1, 1e4, G)          # wildly different mixtures
        m_pairs = ((A @ (kappa * sigma2)) / K_ij
                   + np.array([T[i] + T[j] for i, j in edges]))
        m_random = (W.T.astype(float) @ (kappa * sigma2)) / K_i + T + T_R
        T_hat = retrofit.solve_budgets(
            retrofit.budget_contrasts(m_pairs, m_random, W, kappa), W, kappa)
        truth = np.append(T, T_R); truth = truth - truth.mean()
        assert np.abs(T_hat - truth).max() < 1e-6


def test_report_is_labelled_simulation_only(W):
    rep = retrofit.report(W, n_spots=(100,), seed=0)
    assert 'SIMULATION ONLY' in rep['status']
    assert 'blur' in rep['estimand']          # Lemma D qualification travels
