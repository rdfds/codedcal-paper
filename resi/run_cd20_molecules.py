"""Cross-label statistics on the published RESI molecule positions (one point per
molecule per round), for comparison with the raw-localization analysis.

Usage: python run_cd20_molecules.py CACHE_DIR
"""
import itertools
import json
import sys
from pathlib import Path
import numpy as np
import cd20 as C
from crosslabel.model import bootstrap_fit, fit_dimer
from summarize_cd20 import boot, ci, local_density, merged, pooled_cross

OUT = Path(__file__).resolve().parent / 'results'
EFFICIENCY = 0.5   # labelling efficiency used by Reinhardt et al.


def main(cache):
    rng = np.random.default_rng(2)
    out = {}
    for view in C.VIEWS:
        mols = C.load_resi(cache, view)
        xy = np.concatenate([m['xy'] for m in mols])
        lab = np.concatenate([np.full(len(m['xy']), r) for r, m in enumerate(mols)])
        lp2 = float(np.mean(np.concatenate([m['lp'] for m in mols]) ** 2))
        block = C.blocks(xy)
        N, n = C.block_pair_counts(xy, lab, block=block)
        seed = int(rng.integers(1 << 31))
        obs = C.cross_annuli(xy, lab, seed, fraction=1.0)
        null = np.mean([C.cross_annuli(xy, lab, seed, shift=s, fraction=1.0)
                        for s in C.null_shifts(rng, 12)], 0)
        nd = null[0] / (np.pi * (C.NULL_SMOOTH[1] ** 2 - C.NULL_SMOOTH[0] ** 2))
        pr = obs[1] / null[1]
        rho = local_density(cache, view)
        gfun = lambda a, b: C.g_of_r(pooled_cross(a, b), nd, pr)
        stat = lambda a, b: C.excess_neighbours(gfun(a, b), rho)
        g = gfun(N, n)
        gb = boot(gfun, N, n, rng)
        ex = stat(N, n)
        pairs = {}
        for a, b in itertools.combinations(range(C.K), 2):
            Np, npp = merged(N, n, [[a], [b]])
            pairs[f'R{a+1}R{b+1}'] = [round(stat(Np, npp), 4), ci(boot(stat, Np, npp, rng))]
        se = gb.std(0)
        fit = fit_dimer(C.EDGES, g, se, 2 * lp2)
        bf = bootstrap_fit(C.EDGES, gb[:100], se, 2 * lp2)
        np.savez_compressed(OUT / f'cd20-molecules-{view}.npz', edges=C.EDGES, N=N, n=n, g=g, gboot=gb,
                            null_density=nd, plateau_ratio=pr)
        exb = C.excess_neighbours(gb, rho)
        out[view] = dict(
            treated=C.TREATED[view], molecules=n.sum(0).astype(int).tolist(),
            excess_neighbours_30nm=[round(ex, 4), ci(exb)],
            implied_dimer_subunit_fraction=[round(ex / EFFICIENCY, 3), [round(x / EFFICIENCY, 3) for x in ci(exb)]],
            excess_by_label_pair=pairs,
            dimer_fit=dict(d_nm=[round(fit['params']['d_nm'], 2), ci(bf[:, 1])],
                           sigma_link_nm=[round(fit['params']['sigma_link_nm'], 2), ci(bf[:, 2])],
                           baseline=round(fit['params']['baseline'], 3)),
            g={f'{r:g}nm': round(float(g[int(np.searchsorted((C.EDGES[:-1] + C.EDGES[1:]) / 2, r))]), 3)
               for r in (2.5, 7.5, 12.5, 15.5, 20.5, 30.5, 45.5)})
        print(view, json.dumps(out[view]), flush=True)
    (OUT / 'cd20-molecules.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
