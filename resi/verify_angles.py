"""Check the angle-resolved connected triplet estimator against the true angle
distribution of within-complex triplets, for monomer/dimer mixtures, straight
chains, rings and compact clusters (molecule-level fields, RESI-like precision).

Usage: python verify_angles.py
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import cd20 as C
from simulate import acquire, complexes

SIDE = 10000.0
LINK_SD = 5.0           # nm per molecule, Nup96 fit with the same nanobody label
CASES = [
    dict(name='monomer/dimer', fractions=[0.55, 0.45], shape='chain'),
    dict(name='chains 3-6', fractions=[0.3, 0.2, 0.15, 0.15, 0.1, 0.1], shape='chain'),
    dict(name='hexameric rings', fractions=[0.4, 0.2, 0, 0, 0, 0.4], shape='ring'),
    dict(name='compact clusters 4-6', fractions=[0.4, 0.2, 0, 0.15, 0.15, 0.1], shape='blob'),
]


def oracle_profile(xy, lab, cid, radius=C.TRIPLET_RADIUS, edges=C.ANGLE_EDGES):
    t = cKDTree(xy)
    out = np.zeros(len(edges) - 1)
    for i, js in enumerate(t.query_ball_point(xy, radius)):
        js = np.asarray([j for j in js if j != i and cid[j] == cid[i] and lab[j] != lab[i]], int)
        if len(js) < 2:
            continue
        v = xy[js] - xy[i]; th = np.arctan2(v[:, 1], v[:, 0])
        jj, kk = np.meshgrid(np.arange(len(js)), np.arange(len(js)), indexing='ij')
        ok = lab[js][jj] != lab[js][kk]
        out += np.histogram(np.abs((th[jj] - th[kk] + np.pi) % (2 * np.pi) - np.pi)[ok], edges)[0]
    return out / len(xy)


def run(case, seed):
    rng = np.random.default_rng(seed)
    sub, cid = complexes(rng, SIDE, 150, case['fractions'], 13.5, shape=case['shape'])
    data = acquire(rng, sub, cid, C.K, efficiency=0.5, link_sd=LINK_SD, locs_per_molecule=1, loc_sd=0.0)
    xy = data['site'] + rng.normal(0, 0.5, data['site'].shape)
    lab, cx = data['site_label'], data['site_complex']
    res = C.point_pair_g(xy, lab, rng)
    est, boot = C.triplet_angles(xy, lab, res, rng, B=200)
    truth = oracle_profile(xy, lab, cx)
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
    return dict(case=case['name'], estimate=est.tolist(), lo=lo.tolist(), hi=hi.tolist(), oracle=truth.tolist(),
                covered=float(np.mean((truth >= lo) & (truth <= hi))),
                total_est=float(est.sum()), total_oracle=float(truth.sum()))


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / 'results'
    rows = [run(c, 500 + 10 * k + r) for k, c in enumerate(CASES) for r in range(2)]
    (out / 'verify_angles.json').write_text(json.dumps(dict(angle_edges=C.ANGLE_EDGES.tolist(), rows=rows), indent=1))
    for r in rows:
        print(f"{r['case']:22s} total est {r['total_est']:.4f} oracle {r['total_oracle']:.4f}  bin coverage {r['covered']:.2f}")
        print('   est   ', np.round(r['estimate'], 4).tolist())
        print('   oracle', np.round(r['oracle'], 4).tolist())
