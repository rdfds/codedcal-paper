"""Sensitivity of the main results to analysis choices.

CD20 (published RESI molecule positions):
  - excess labelled neighbours: integration radius x baseline window
  - connected triplets m2: radius x local-background annulus
  - interval width: spatial block size
NPC (Nup96):
  - pore-level efficiency: molecule-to-pore assignment radius
  - corner-level efficiency: radius, linkage spread and corner geometry

Usage: python robustness.py CACHE_DIR
"""
import itertools
import json
import sys
from pathlib import Path
import numpy as np
import cd20 as C
import npc as P
from npc_validation import load_npc

OUT = Path(__file__).resolve().parent / 'results'
B = 200
RMAX = [20.0, 25.0, 30.0, 40.0]
BASES = [(50.0, 90.0), (60.0, 100.0), (70.0, 100.0)]
TRIPLET_R = [20.0, 25.0, 30.0, 40.0]
ANNULI = [(80.0, 200.0), (100.0, 250.0), (150.0, 300.0)]
BLOCKS = [1000.0, 2000.0, 4000.0]


def ci(x):
    return np.round(np.percentile(x, [2.5, 97.5]), 4).tolist()


def cd20_view(cache, view, rng):
    mols = C.load_resi(cache, view)
    xy = np.concatenate([m['xy'] for m in mols])
    lab = np.concatenate([np.full(len(m['xy']), r) for r, m in enumerate(mols)])
    rho = C.plateau_density(xy)
    out = dict(excess={}, m2={}, block_ci_width={})
    res = C.point_pair_g(xy, lab, rng)
    nd, pr = res['null_density'], res['plateau_ratio']
    for rmax, base in itertools.product(RMAX, BASES):
        fn = lambda N, n: C.excess_neighbours(C.g_of_r(C.cross_g(N, n)[1], nd, pr), rho, rmax=rmax, base=base)
        out['excess'][f'r{rmax:g}_b{base[0]:g}-{base[1]:g}'] = [round(float(fn(res['N'], res['n'])), 4),
                                                                ci(C.bootstrap(fn, res['N'], res['n'], rng, B))]
    for radius, ann in itertools.product(TRIPLET_R, ANNULI):
        m2, m2ci, _ = C.triplet_m2(xy, lab, res, rng, B=B, radius=radius, annulus=ann)
        out['m2'][f'r{radius:g}_a{ann[0]:g}-{ann[1]:g}'] = [round(m2, 4), np.round(m2ci, 4).tolist()]
    for bs in BLOCKS:
        rb = C.point_pair_g(xy, lab, rng, block_size=bs)
        val, (lo, hi) = C.excess_with_ci(rb, rho, rng, B)
        out['block_ci_width'][f'{bs:g}'] = dict(excess=round(val, 4), width=round(hi - lo, 4),
                                                 blocks=int(rb['block'].max() + 1))
    return out


def npc_sensitivity(cache, rng):
    out = dict(assignment={}, corner={})
    for assign in (70.0, 90.0, 110.0):
        xy, lab, pid, info = load_npc(cache, assign)
        pi = np.bincount(lab, minlength=P.K) / len(lab)
        st = P.stoichiometry(P.pore_counts(lab, pid), pi)
        out['assignment'][f'{assign:g}'] = {k: round(float(v), 4) for k, v in st.items()} | dict(assigned=info['assigned'])
    xy, lab, pid, _ = load_npc(cache)
    pi = np.bincount(lab, minlength=P.K) / len(lab)
    defaults = (P.RING_RADIUS, P.PAIR_SPACING, P.RING_OFFSET)
    stats = {}
    for radius in (20.0, 22.0, 24.0, 26.0, 28.0):
        c, t = P.corner_statistics(xy, lab, pid, pi, radius)
        stats[radius] = (c.mean(), t.mean())
    grid = list(itertools.product((20.0, 24.0, 28.0), (4.0, 5.0, 6.0), (49.0, 51.0, 53.0), (10.0, 12.0, 14.0),
                                  (0.0, 3.0, 6.0)))
    for radius, sd, ring, spacing, offset in grid:
        P.RING_RADIUS, P.PAIR_SPACING, P.RING_OFFSET = ring, spacing, offset
        a1, a2 = P.geometry_factors(rng, radius, link_sd=sd, n_pores=150)
        nex, m2 = stats[radius]
        out['corner'][f'R{radius:g}_sd{sd:g}_ring{ring:g}_pair{spacing:g}_off{offset:g}'] = dict(
            p_pairs=round(nex / a1, 4), p_triplets=round(float(np.sqrt(m2 / a2)), 4),
            p_ratio=round((m2 / nex) * (a1 / a2), 4))
    P.RING_RADIUS, P.PAIR_SPACING, P.RING_OFFSET = defaults
    vals = np.array([[v['p_pairs'], v['p_triplets'], v['p_ratio']] for v in out['corner'].values()])
    out['corner_range'] = dict(p_pairs=[float(vals[:, 0].min()), float(vals[:, 0].max())],
                               p_triplets=[float(vals[:, 1].min()), float(vals[:, 1].max())],
                               p_ratio=[float(vals[:, 2].min()), float(vals[:, 2].max())])
    return out


def main(cache):
    rng = np.random.default_rng(7)
    out = dict(cd20={}, npc=None)
    for view in C.VIEWS:
        out['cd20'][view] = cd20_view(cache, view, rng)
        ex = [v[0] for v in out['cd20'][view]['excess'].values()]
        m2 = [v[0] for v in out['cd20'][view]['m2'].values()]
        print(view, 'excess range', round(min(ex), 3), round(max(ex), 3), '| m2 range', round(min(m2), 3),
              round(max(m2), 3), '| block widths',
              {k: v['width'] for k, v in out['cd20'][view]['block_ci_width'].items()}, flush=True)
    out['npc'] = npc_sensitivity(cache, rng)
    print('npc assignment', out['npc']['assignment'])
    print('npc corner range', out['npc']['corner_range'])
    (OUT / 'robustness.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
