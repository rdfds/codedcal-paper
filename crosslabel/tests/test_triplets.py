import numpy as np

from crosslabel.correlation import point_pair_g, triplet_m2
from crosslabel.labels import connected_triplet_angles, connected_triplets


def test_monomer_dimer_has_no_connected_triplets(molecule_field):
    xy, lab, _, rng = molecule_field([0.55, 0.45], seed=2)
    m2, ci, _ = triplet_m2(xy, lab, point_pair_g(xy, lab, rng), rng, B=100)
    assert abs(m2) < 0.03 and ci[0] < 0.0 < ci[1] + 0.02


def test_trimers_are_detected(molecule_field):
    xy, lab, _, rng = molecule_field([0.4, 0.3, 0.3], seed=3)
    m2, ci, _ = triplet_m2(xy, lab, point_pair_g(xy, lab, rng), rng, B=100)
    assert ci[0] > 0.02


def test_angle_profile_sums_to_triplet_statistic(molecule_field):
    xy, lab, _, _ = molecule_field([0.4, 0.3, 0.3], seed=4)
    edges = np.linspace(0, np.pi, 7)
    zero = np.zeros((4, 4)); zero_a = np.zeros((4, 4, len(edges) - 1))
    q = connected_triplets(xy, lab, 4, 30.0, (100.0, 250.0), zero)
    qa = connected_triplet_angles(xy, lab, 4, 30.0, (100.0, 250.0), zero_a, edges)
    assert np.allclose(qa.sum(1), q)
