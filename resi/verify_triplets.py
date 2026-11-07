"""Check the connected triplet estimator on simulated molecule-level fields.

Each field gives one point per labelled molecule (RESI-like precision). The
estimator m2 = mean(q) / (2 P3) should equal the mean number of unordered pairs
of labelled complex partners within R per labelled molecule: zero for any
monomer/dimer mixture (with or without large-scale density variation) and
positive for larger complexes.
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import cd20 as C
from crosslabel.simulate import acquire, complexes, smooth_field

SIDE = 10000.0
RADIUS = C.TRIPLET_RADIUS
CASES = [
    dict(name='monomer/dimer', fractions=[0.55, 0.45], hetero=False),
    dict(name='monomer/dimer, heterogeneous', fractions=[0.55, 0.45], hetero=True),
    dict(name='monomer only, heterogeneous', fractions=[1.0], hetero=True),
    dict(name='with trimers', fractions=[0.4, 0.3, 0.3], hetero=False),
    dict(name='chains up to 5, heterogeneous', fractions=[0.3, 0.2, 0.2, 0.15, 0.15], hetero=True),
]


def oracle_m2(xy, cid, radius=RADIUS):
    t = cKDTree(xy)
    nb = t.query_ball_point(xy, radius)
    c = np.array([sum(cid[j] == cid[i] for j in js if j != i) for i, js in enumerate(nb)])
    return float((c * (c - 1) / 2).mean())


def run(case, seed):
    rng = np.random.default_rng(seed)
    field = smooth_field(rng, SIDE) if case['hetero'] else None
    sub, cid = complexes(rng, SIDE, 150, case['fractions'], 13.5, intensity=field)
    data = acquire(rng, sub, cid, C.K, efficiency=0.5, link_sd=4.0, locs_per_molecule=1, loc_sd=0.0)
    xy = data['site'] + rng.normal(0, 0.5, data['site'].shape)
    lab, cx = data['site_label'], data['site_complex']
    m2, ci, p3 = C.triplet_m2(xy, lab, C.point_pair_g(xy, lab, rng), rng, B=200)
    truth = oracle_m2(xy, cx)
    return dict(case=case['name'], molecules=len(xy), m2=m2, ci=ci, oracle=truth,
                covered=bool(ci[0] <= truth <= ci[1]), p3=p3)


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / 'results'
    rows = [run(c, 300 + 10 * k + r) for k, c in enumerate(CASES) for r in range(3)]
    (out / 'verify_triplets.json').write_text(json.dumps(rows, indent=1))
    for r in rows:
        print(f"{r['case']:32s} n={r['molecules']:6d}  m2={r['m2']:.4f}  "
              f"CI=[{r['ci'][0]:.4f},{r['ci'][1]:.4f}]  oracle={r['oracle']:.4f}  covered={r['covered']}")
    print('coverage', np.mean([r['covered'] for r in rows]))
