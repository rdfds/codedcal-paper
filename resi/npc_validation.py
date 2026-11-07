"""Nup96 NPC validation of the label statistics (Reinhardt et al. 2023, Fig. 2 data).

1. Simulation: pores with known efficiency p and the observed label frequencies,
   with and without same-label cluster merging; every route should recover p.
2. Linkage spread: chosen to match the observed same-pore cross-label
   nearest-neighbour distance distribution (corner and adjacent-corner peaks).
3. Data: p from molecule counts, cross-label pairs and three-distinct-label
   triples per pore (n = 32), (n, p) jointly from pairs and triples, and corner-
   level estimates within the inter-corner gap.

Usage: python npc_validation.py CACHE_DIR
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import npc as P

OUT = Path(__file__).resolve().parent / 'results'
RADIUS = 24.0            # nm, the gap between corner (~12 nm) and adjacent-corner (~41 nm) peaks
ASSIGN_NM = 90.0         # molecule-to-pore-centre assignment radius
NN_EDGES = np.arange(0, 64, 4.0)
B = 400


def load_npc(cache):
    xy, lab, cen = [], [], {}
    for r in range(1, 5):
        d = np.load(Path(cache) / f'npc-resi-R{r}.npz')
        xy.append(np.column_stack([d['x'], d['y']])); lab.append(np.full(len(d['x']), r - 1))
        loc = np.load(Path(cache) / f'npc-locs-R{r}.npz')
        for g in np.unique(loc['group']):
            m = loc['group'] == g
            s = cen.setdefault(int(g), [0.0, 0.0, 0])
            s[0] += float(loc['x'][m].sum()); s[1] += float(loc['y'][m].sum()); s[2] += int(m.sum())
    groups = np.array(sorted(cen))
    centres = np.array([[cen[g][0] / cen[g][2], cen[g][1] / cen[g][2]] for g in groups])
    xy, lab = np.concatenate(xy), np.concatenate(lab)
    d, idx = cKDTree(centres).query(xy, distance_upper_bound=ASSIGN_NM)
    ok = np.isfinite(d)
    return xy[ok], lab[ok], idx[ok], dict(molecules=int(len(xy)), assigned=int(ok.sum()), pores=int(len(groups)))


def nn_hist(xy, lab, pid):
    """Same-pore nearest cross-label neighbour distance histogram (density)."""
    out = []
    for a in range(P.K):
        for b in range(P.K):
            if a == b:
                continue
            A = np.flatnonzero(lab == a); Bi = np.flatnonzero(lab == b)
            if not len(A) or not len(Bi):
                continue
            t = cKDTree(xy[Bi])
            d, j = t.query(xy[A], k=8, distance_upper_bound=NN_EDGES[-1])
            for row_d, row_j, i in zip(d, j, A):
                okm = np.isfinite(row_d)
                same = okm.copy()
                same[okm] = pid[Bi[row_j[okm]]] == pid[i]
                if same.any():
                    out.append(row_d[same][0])
    h, _ = np.histogram(out, NN_EDGES)
    return h / max(h.sum(), 1)


def estimate(xy, lab, pid, pi, a1, a2, rng):
    N = P.pore_counts(lab, pid)
    st = P.stoichiometry(N, pi)
    c, t = P.corner_statistics(xy, lab, pid, pi, RADIUS)
    nex, m2 = c.mean(), t.mean()
    out = dict(st, n_ex=nex, m2=m2, p_corner_pairs=nex / a1, p_corner_triplets=np.sqrt(max(m2, 0) / a2),
               p_corner_ratio=(m2 / nex) * (a1 / a2))
    # bootstrap over pores
    u, inv = np.unique(pid, return_inverse=True)
    rows = []
    cs = np.bincount(inv, weights=c, minlength=len(u)); ts = np.bincount(inv, weights=t, minlength=len(u))
    cnt = np.bincount(inv, minlength=len(u)).astype(float)
    for _ in range(B):
        k = rng.integers(0, len(u), len(u))
        stb = P.stoichiometry(N[k], pi)
        nb_, mb_ = cs[k].sum() / cnt[k].sum(), ts[k].sum() / cnt[k].sum()
        rows.append([stb['p_counts'], stb['p_pairs'], stb['p_triples'], stb['n_joint'], stb['p_joint'],
                     nb_ / a1, np.sqrt(max(mb_, 0) / a2), (mb_ / nb_) * (a1 / a2)])
    keys = ['p_counts', 'p_pairs', 'p_triples', 'n_joint', 'p_joint', 'p_corner_pairs',
            'p_corner_triplets', 'p_corner_ratio']
    rows = np.array(rows)
    out['ci'] = {k: np.nanpercentile(rows[:, i], [2.5, 97.5]).tolist() for i, k in enumerate(keys)}
    return {k: (float(v) if not isinstance(v, dict) else v) for k, v in out.items()}


def main(cache):
    rng = np.random.default_rng(6)
    xy, lab, pid, info = load_npc(cache)
    pi = np.bincount(lab, minlength=P.K) / len(lab)
    obs_h = nn_hist(xy, lab, pid)

    # linkage spread matching the observed corner structure
    fits = {}
    for sd in (2.0, 3.0, 4.0, 5.0, 6.0, 7.0):
        sx, sl, sp = P.simulate(rng, 600, 0.4, pi, link_sd=sd)
        fits[sd] = float(((P_nn := nn_hist(sx, sl, sp)) - obs_h) @ (P_nn - obs_h))
    link_sd = min(fits, key=fits.get)
    a1, a2 = P.geometry_factors(rng, RADIUS, link_sd=link_sd)
    sens = {f'{sd:g}': P.geometry_factors(rng, RADIUS, link_sd=sd, n_pores=200) for sd in (link_sd - 1, link_sd + 1)}

    # simulation recovery
    sims = []
    for p_true in (0.3, 0.4, 0.5):
        for merge in (0.0, 8.0):
            sx, sl, sp = P.simulate(rng, 1200, p_true, pi, link_sd=link_sd, merge_nm=merge)
            spi = np.bincount(sl, minlength=P.K) / len(sl)
            e = estimate(sx, sl, sp, spi, a1, a2, rng)
            sims.append(dict(p_true=p_true, merge_nm=merge, **{k: v for k, v in e.items() if k != 'ci'}, ci=e['ci']))
            print('sim p=%.2f merge=%g' % (p_true, merge),
                  {k: round(e[k], 3) for k in ('p_counts', 'p_pairs', 'p_triples', 'n_joint',
                                              'p_corner_pairs', 'p_corner_triplets', 'p_corner_ratio')}, flush=True)

    data = estimate(xy, lab, pid, pi, a1, a2, rng)
    print('data', {k: round(v, 3) for k, v in data.items() if k != 'ci'})
    print('data CIs', {k: np.round(v, 3).tolist() for k, v in data['ci'].items()})
    out = dict(info=info, label_fractions=pi.round(4).tolist(), radius_nm=RADIUS,
               link_sd_fit=link_sd, link_sd_misfit=fits, a1=a1, a2=a2, geometry_sensitivity=sens,
               observed_nn_hist=obs_h.tolist(), nn_edges=NN_EDGES.tolist(),
               simulations=sims, data=data)
    (OUT / 'npc-validation.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
