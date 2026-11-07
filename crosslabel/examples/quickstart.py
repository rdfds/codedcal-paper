"""Quickstart: simulate a stochastically labelled field and measure its organization.

Two fields with identical pair statistics are generated (dimers vs equilateral
trimers holding half as many subunits); pair statistics agree, the connected
triplet statistic separates them, and the angle profile shows their shapes.

Run: python examples/quickstart.py
"""
import numpy as np
import crosslabel as cl
from crosslabel.simulate import acquire, complexes

K = 4                    # labels (e.g. DNA-PAINT docking sequences / rounds)
EFFICIENCY = 0.5         # fraction of subunits labelled and detected


def field(fractions, shape, seed):
    rng = np.random.default_rng(seed)
    subunits, complex_id = complexes(rng, side=10000.0, density=100.0, fractions=fractions,
                                     spacing=13.5, shape=shape)
    data = acquire(rng, subunits, complex_id, K, efficiency=EFFICIENCY, link_sd=5.0,
                   locs_per_molecule=1, loc_sd=0.0)
    xy = data['site'] + rng.normal(0, 0.5, data['site'].shape)   # one point per molecule (nm)
    return xy, data['site_label'], rng


for name, fractions, shape in (('dimers', [1 - 0.6 / 1.4, 0.6 / 1.4], 'chain'),
                               ('triangular trimers', [0.875, 0.0, 0.125], 'ring')):
    xy, labels, rng = field(fractions, shape, seed=1)
    res = cl.cross_label_g(xy, labels, rng)                     # pair counts, null, g(r)
    rho = cl.plateau_density(xy)
    excess, excess_ci = cl.excess_with_ci(res, rho, rng)        # labelled partners within 30 nm
    m2, m2_ci, _ = cl.connected_triplets_m2(xy, labels, res, rng)
    angles, _ = cl.triplet_angle_profile(xy, labels, res, rng)
    print(f'{name:20s} excess neighbours {excess:.2f} [{excess_ci[0]:.2f}, {excess_ci[1]:.2f}]   '
          f'm2 {m2:.3f} [{m2_ci[0]:.3f}, {m2_ci[1]:.3f}]')
    if m2_ci[0] > 0:   # shape is only defined when connected triplets are present
        print(f'{"":20s} angle profile (0-180 deg, 15-deg bins): {np.round(angles / angles.sum(), 2).tolist()}')
    else:
        print(f'{"":20s} no connected triplets: no higher-order assemblies to shape')

# planning: probability that a trimer's three subunits carry distinct labels, by number of labels
print('three distinct labels:', {k: round(cl.design.uniform_distinct(3, k), 3) for k in (3, 4, 6, 8, 12)})
