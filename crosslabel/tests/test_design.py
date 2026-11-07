import numpy as np

from crosslabel.design import distinct_probability, fixed_budget, uniform_distinct


def test_distinct_probabilities():
    assert np.isclose(uniform_distinct(2, 4), 0.75)
    assert np.isclose(uniform_distinct(3, 4), 0.375)
    for K in (3, 4, 8):
        for m in (2, 3):
            assert np.isclose(distinct_probability(m, np.ones(K)), uniform_distinct(m, K))
    assert uniform_distinct(3, 2) == 0.0


def test_fixed_budget_reference_is_one():
    out = fixed_budget(3, [4, 8], lambda phi: min(1.0, phi))
    assert np.isclose(out[4], 1.0)
