"""Per-view cross-label analysis of the CD20 acquisitions.

Usage: python run_cd20.py CACHE_DIR [view ...]
Writes results/cd20-<view>.npz and results/cd20-summary.json.
"""
import json
import sys
from pathlib import Path
import numpy as np
import cd20 as C

OUT = Path(__file__).resolve().parent / 'results'
SHORT = 40.0   # nm, window for same-molecule pairs


def analyse(cache, view, seed=0):
    rng = np.random.default_rng(seed)
    d = C.load(cache, view)
    xy, lab = d['xy'], d['label']
    block = C.blocks(xy)
    N, n = C.block_pair_counts(xy, lab, block=block)
    null_density, plateau_ratio = C.null_normalization(xy, lab, rng)
    G, pooled = C.cross_g(N, n)
    g = C.g_of_r(pooled, null_density, plateau_ratio)

    def g_fn(Nb, nb):
        return C.g_of_r(C.cross_g(Nb, nb)[1], null_density, plateau_ratio)
    gboot = C.bootstrap(g_fn, N, n, rng)

    # label exchangeability: each label pair's g with the common normalization
    pair_g = np.array([[C.g_of_r(G[a, b], null_density, plateau_ratio) if a != b
                        else np.full(len(g), np.nan) for b in range(C.K)] for a in range(C.K)])

    S, _ = C.decomposition(N, n)
    centres = (C.EDGES[:-1] + C.EDGES[1:]) / 2
    short = centres < SHORT
    ns = n.sum(0)
    s_tot = S[:, short].sum(1)
    poisson_molecules = ns ** 2 / s_tot
    resi = C.load_resi(cache, view)
    resi_counts = np.array([len(r['xy']) for r in resi])
    offsets, offset_se, close_pairs = C.registration_offsets(xy, lab)

    np.savez_compressed(OUT / f'cd20-{view}.npz', edges=C.EDGES, N=N, n=n,
                        null_density=null_density, plateau_ratio=plateau_ratio,
                        g=g, gboot=gboot, pair_g=pair_g, S=S, offsets=offsets)
    lo, hi = np.percentile(gboot, [2.5, 97.5], axis=0)
    pick = lambda r: int(np.searchsorted(centres, r))
    return dict(
        view=view, treated=C.TREATED[view], localizations=ns.astype(int).tolist(),
        mean_lp_nm=[float(d['lp'][lab == a].mean()) for a in range(C.K)],
        blocks=int(N.shape[0]),
        g={f'{r:g}nm': [float(g[pick(r)]), float(lo[pick(r)]), float(hi[pick(r)])]
           for r in (2.5, 7.5, 12.5, 17.5, 25.5, 40.5, 57.5, 92.5)},
        poisson_molecules=poisson_molecules.round().astype(int).tolist(),
        resi_molecules=resi_counts.tolist(),
        locs_per_molecule_poisson=(s_tot / ns).round(1).tolist(),
        registration_offset_nm={f'R{a+1}-R{b+1}': [np.round(offsets[a, b], 2).tolist(),
                                                    np.round(offset_se[a, b], 2).tolist()]
                                for a in range(C.K) for b in range(a + 1, C.K)},
        plateau_ratio=float(plateau_ratio),
    )


if __name__ == '__main__':
    cache = Path(sys.argv[1])
    views = sys.argv[2:] or C.VIEWS
    OUT.mkdir(exist_ok=True)
    path = OUT / 'cd20-summary.json'
    summary = json.loads(path.read_text()) if path.exists() else {}
    for v in views:
        summary[v] = analyse(cache, v)
        print(json.dumps(summary[v], indent=1), flush=True)
        path.write_text(json.dumps(summary, indent=1))
