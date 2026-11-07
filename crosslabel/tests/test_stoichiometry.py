import numpy as np

from crosslabel.stoichiometry import group_counts, stoichiometry


def test_efficiency_from_pairs_and_triples():
    rng = np.random.default_rng(5)
    n, p, pi = 32, 0.4, np.array([0.14, 0.25, 0.42, 0.19])
    groups, labels = [], []
    for g in range(3000):
        k = rng.binomial(n, p)
        labels.append(rng.choice(4, k, p=pi)); groups.append(np.full(k, g))
    lab, grp = np.concatenate(labels), np.concatenate(groups)
    est = stoichiometry(group_counts(lab, grp, 4), np.bincount(lab) / len(lab), n)
    for key in ('p_counts', 'p_pairs', 'p_triples'):
        assert abs(est[key] - p) < 0.01, key
    assert 24 < est['n_joint'] < 42
