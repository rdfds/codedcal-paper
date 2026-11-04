"""Exact finite-law checks for clarified theoretical examples; no sampled study."""
from collections import Counter
from fractions import Fraction as F
from itertools import product


def test_fourth_order_all_source_partitions():
    coefficient = {1:F(1), 2:F(-1,3), 3:F(1,3), 4:F(-1)}
    partitions = ((0,0,0,0), (0,0,0,1), (0,0,1,1),
                  (0,0,1,2), (0,1,2,3))
    distributions = []
    for partition in partitions:
        sources = len(set(partition))
        counts = Counter(len({codes[i] for i in partition})
                         for codes in product(range(4), repeat=sources))
        law = {d:F(n,4**sources) for d,n in counts.items()}
        assert sum(coefficient[d]*p for d,p in law.items()) == (1 if sources==1 else 0)
        distributions.append(law)
    assert distributions[1] == distributions[2] == {1:F(1,4),2:F(3,4)}


def test_count_toy_asymptotic_factor_is_not_a_finite_ratio_variance():
    values = []
    for code in product(range(4), repeat=3):
        count = [code.count(c) for c in range(4)]
        values.append(F(16*sum(n**3 for n in count)-36*sum(n*n for n in count)+54,6))
    law = {w:F(n,64) for w,n in Counter(values).items()}
    assert law == {F(27):F(1,16),F(3):F(9,16),F(-1):F(6,16)}
    mean = sum(w*p for w,p in law.items())
    second = sum(w*w*p for w,p in law.items())
    assert mean==3 and second==51
    assert second/mean**2==F(17,3)
    assert sum((w-mean)**4*p for w,p in law.items())==20832
    # Four iid groups have positive probability of a zero pooled denominator.
    assert 3 + (-1)*3 == 0
    assert 4*law[F(3)]*law[F(-1)]**3 > 0
