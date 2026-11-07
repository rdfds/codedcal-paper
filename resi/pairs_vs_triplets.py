"""Two oligomer models that pair statistics cannot tell apart.

Model D: monomers and dimers (spacing d), 60% of subunits in dimers.
Model T: monomers and equilateral trimers (side d), 30% of subunits in trimers.

Every partner sits at distance d in both models, and each subunit has on average
0.6 partners in both, so the partner-distance distribution and pair excess are
identical: the cross-label pair correlation g(r), excess neighbours and any
other second-order statistic coincide in expectation. Only the trimers have
pairs of partners, so the connected triplet statistic is zero for D and
positive for T.

Usage: python pairs_vs_triplets.py
"""
import json
import sys
from pathlib import Path
import numpy as np
import cd20 as C
from simulate import acquire, complexes

SIDE = 10000.0
DENSITY = 100.0
SPACING = 13.5
LINK_SD = 5.0
EFFICIENCY = 0.5
MODELS = {
    'dimers (60% of subunits)': dict(fractions=[1 - 0.6 / 1.4, 0.6 / 1.4], shape='chain'),
    'triangular trimers (30%)': dict(fractions=[0.875, 0.0, 0.125], shape='ring'),
}


def run(name, seed):
    rng = np.random.default_rng(seed)
    m = MODELS[name]
    sub, cid = complexes(rng, SIDE, DENSITY, m['fractions'], SPACING, shape=m['shape'])
    d = acquire(rng, sub, cid, C.K, efficiency=EFFICIENCY, link_sd=LINK_SD, locs_per_molecule=1, loc_sd=0.0)
    xy = d['site'] + rng.normal(0, 0.5, d['site'].shape)
    lab = d['site_label']
    res = C.point_pair_g(xy, lab, rng)
    rho = C.plateau_density(xy)
    ex, ex_ci = C.excess_with_ci(res, rho, rng, B=200)
    m2, m2_ci, _ = C.triplet_m2(xy, lab, res, rng, B=200)
    return dict(model=name, g=res['g'].tolist(), excess=ex, excess_ci=ex_ci, m2=m2, m2_ci=m2_ci,
                expected_m2=EFFICIENCY ** 2 * (0.3 if 'trimer' in name else 0.0))


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / 'results'
    rows = [run(name, 700 + 10 * k + r) for k, name in enumerate(MODELS) for r in range(4)]
    (out / 'pairs-vs-triplets.json').write_text(json.dumps(dict(edges=C.EDGES.tolist(), rows=rows), indent=1))
    for name in MODELS:
        rr = [r for r in rows if r['model'] == name]
        print(f"{name:26s} excess {np.mean([r['excess'] for r in rr]):.3f} (runs {[round(r['excess'], 3) for r in rr]})"
              f" | m2 {np.mean([r['m2'] for r in rr]):.4f} (runs {[round(r['m2'], 4) for r in rr]},"
              f" expected {rr[0]['expected_m2']:.3f})")
    gd = np.mean([r['g'] for r in rows if 'dimer' in r['model']], 0)
    gt = np.mean([r['g'] for r in rows if 'trimer' in r['model']], 0)
    c = (C.EDGES[:-1] + C.EDGES[1:]) / 2
    for rr in (2.5, 7.5, 12.5, 17.5, 25.5, 40.5):
        k = int(np.searchsorted(c, rr))
        print(f'  g({c[k]:.1f} nm): dimers {gd[k]:.3f}  trimers {gt[k]:.3f}')
