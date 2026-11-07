import numpy as np

from crosslabel.labels import decompose, label_moments, pair_counts, split_pairs, star_triplets
from crosslabel.oracle import exact_cross_variance, oracle_pairs, oracle_star
from crosslabel.simulate import acquire, complexes, relabel

EDGES = np.arange(0, 31, 10.0)


def test_label_moments_uniform():
    m = label_moments(np.ones(4))
    assert np.isclose(m['cross'], 0.75) and np.isclose(m['distinct3'], 0.375)


def test_exact_identities_under_label_redraws():
    rng = np.random.default_rng(1)
    K, pi = 4, np.array([0.4, 0.3, 0.2, 0.1])
    sub, cid = complexes(rng, 1200, 120, [0.5, 0.3, 0.2], 13.5)
    data = acquire(rng, sub, cid, K, efficiency=0.6, pi=pi, link_sd=3.0, locs_per_molecule=8, loc_sd=4.0)
    xy, mol = data['xy'], data['mol']
    S, W, pairs = oracle_pairs(xy, mol, EDGES)
    W3 = oracle_star(xy, mol, 20.0)
    m = label_moments(pi)
    draws = [relabel(rng, data, K, pi) for _ in range(300)]
    cross = np.array([split_pairs(pair_counts(xy, lab, EDGES, K))[1] for lab in draws])
    Shat = np.array([decompose(pair_counts(xy, lab, EDGES, K), pi)['S'] for lab in draws])
    star = np.array([star_triplets(xy, lab, 20.0, K) for lab in draws])
    z = lambda x, e: (x.mean(0) - e) / (x.std(0, ddof=1) / np.sqrt(len(x)))
    assert np.all(np.abs(z(cross, m['cross'] * W)) < 4.5)
    assert np.all(np.abs(z(Shat, S)) < 4.5)
    assert abs(z(star, m['distinct3'] * W3)) < 4.5
    ratio = cross.var(0, ddof=1) / exact_cross_variance(pairs, pi, len(EDGES) - 1, int(mol.max()) + 1)
    assert np.all((ratio > 0.7) & (ratio < 1.4))
