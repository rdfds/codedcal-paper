"""Design checks on the CD20 data.

1. Round length: raw localizations are reassigned to the published RESI molecules
   (same round, nearest centre within 9 nm). Keeping only the first fraction phi
   of each round's frames, a molecule is detected when it retains at least
   MIN_LOCS localizations; its position is the precision-weighted mean of the
   retained localizations. Pair and triplet statistics are recomputed on the
   detected molecules. If detection is independent of organization, excess
   neighbours scale with p_det and m2 with p_det^2.
2. Number of labels: the four rounds are merged into three or two labels on the
   molecule positions; interval widths are compared with the label factors
   1 - s2 (pairs) and the three-distinct-label probability (triplets).
3. Label exchangeability: the spread of excess neighbours across the six label
   pairs is compared with its distribution under random permutation of labels
   over fixed molecule positions.

Usage: python design_cd20.py CACHE_DIR
"""
import itertools
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import cd20 as C
from crosslabel.design import fixed_budget, interpolated
from crosslabel.labels import label_moments

OUT = Path(__file__).resolve().parent / 'results'
PHIS = [1.0, 0.5, 0.25, 0.125, 0.0625]
MIN_LOCS = 16          # smallest localization count in the published RESI tables
ASSIGN_NM = 9.0        # RESI clustering radius
PERMUTATIONS = 500


def molecules_from_rounds(cache, view):
    """Reassign raw localizations to RESI molecules; returns per-molecule arrays."""
    d = C.load(cache, view)
    resi = C.load_resi(cache, view)
    out = []
    for r, m in enumerate(resi):
        sel = d['label'] == r
        xy, fr, lp = d['xy'][sel], d['frame'][sel], d['lp'][sel]
        dist, idx = cKDTree(m['xy']).query(xy, distance_upper_bound=ASSIGN_NM)
        ok = np.isfinite(dist)
        out.append(dict(xy=xy[ok], frame=fr[ok], w=1 / lp[ok] ** 2, mol=idx[ok],
                        n_mol=len(m['xy']), frames=int(d['frame'][sel].max()) + 1))
    return out


def detect(rounds, phi, min_locs=MIN_LOCS):
    xy, lab = [], []
    for r, d in enumerate(rounds):
        keep = d['frame'] < phi * d['frames']
        mol, w = d['mol'][keep], d['w'][keep]
        n = np.bincount(mol, minlength=d['n_mol'])
        sw = np.bincount(mol, weights=w, minlength=d['n_mol'])
        pos = np.column_stack([np.bincount(mol, weights=w * d['xy'][keep, c], minlength=d['n_mol'])
                               for c in range(2)]) / np.maximum(sw, 1e-12)[:, None]
        hit = n >= min_locs
        xy.append(pos[hit]); lab.append(np.full(hit.sum(), r))
    return np.concatenate(xy), np.concatenate(lab)


def statistics(xy, lab, rng, triplets=True):
    res = C.point_pair_g(xy, lab, rng)
    rho = C.plateau_density(xy)
    ex, ex_ci = C.excess_with_ci(res, rho, rng)
    out = dict(molecules=int(len(xy)), excess=ex, excess_ci=ex_ci)
    if triplets and lab.max() >= 2:
        m2, m2_ci, p3 = C.triplet_m2(xy, lab, res, rng)
        out.update(m2=m2, m2_ci=m2_ci, p3=p3)
    return out


def pair_excess_spread(pairs_d, li, lj, lab, rho, edges=C.EDGES):
    """SD across label pairs of the self-normalized excess neighbours."""
    c = (edges[:-1] + edges[1:]) / 2
    ann = np.pi * np.diff(edges ** 2)
    base = (c >= 60) & (c < 100)
    inner = edges[1:] <= 30
    b = np.digitize(pairs_d, edges) - 1
    a, bb = np.minimum(lab[li], lab[lj]), np.maximum(lab[li], lab[lj])
    vals = []
    for p, q in itertools.combinations(range(C.K), 2):
        s = (a == p) & (bb == q)
        h = np.bincount(b[s], minlength=len(ann))[:len(ann)] / ann
        g = h / (h[base] * ann[base]).sum() * ann[base].sum()
        vals.append(rho * ((g[inner] - 1) * ann[inner]).sum())
    return float(np.std(vals)), vals


def main(cache):
    rng = np.random.default_rng(4)
    out = {}
    for view in C.VIEWS:
        rounds = molecules_from_rounds(cache, view)
        n_resi = sum(d['n_mol'] for d in rounds)
        thinning = {}
        for phi in PHIS:
            xy, lab = detect(rounds, phi)
            s = statistics(xy, lab, rng)
            s['p_det'] = len(xy) / n_resi
            thinning[str(phi)] = s
        # label merging on the published molecule positions
        mols = C.load_resi(cache, view)
        mxy = np.concatenate([m['xy'] for m in mols])
        mlab = np.concatenate([np.full(len(m['xy']), r) for r, m in enumerate(mols)])
        merging = {}
        for k in (4, 3, 2):
            rows = []
            for groups in partitions(k):
                relab = np.zeros_like(mlab)
                for gi, g in enumerate(groups):
                    relab[np.isin(mlab, g)] = gi
                s = statistics(mxy, relab, rng, triplets=k >= 3)
                pi = np.bincount(relab, minlength=k) / len(relab)
                m = label_moments(pi)
                s.update(groups=groups, cross=m['cross'], distinct3=m['distinct3'])
                rows.append(s)
            merging[f'K={k}'] = rows
        # exchangeability permutation test on molecule positions
        t = cKDTree(mxy)
        li, lj = t.query_pairs(C.EDGES[-1], output_type='ndarray').T
        dd = np.hypot(*(mxy[li] - mxy[lj]).T)
        rho = C.plateau_density(mxy)
        obs, vals = pair_excess_spread(dd, li, lj, mlab, rho)
        null = np.array([pair_excess_spread(dd, li, lj, rng.permutation(mlab), rho)[0]
                         for _ in range(PERMUTATIONS)])
        out[view] = dict(treated=C.TREATED[view], resi_molecules=n_resi, thinning=thinning,
                         merging=merging,
                         exchangeability=dict(pair_excess=vals, sd=obs,
                                              null_sd_median=float(np.median(null)),
                                              p_value=float((1 + (null >= obs).sum()) / (1 + len(null)))))
        th = thinning
        print(view, 'p_det', [round(th[str(p)]['p_det'], 3) for p in PHIS],
              '| excess', [round(th[str(p)]['excess'], 3) for p in PHIS],
              '| m2', [round(th[str(p)].get('m2', np.nan), 3) for p in PHIS],
              '| exch p', round(out[view]['exchangeability']['p_value'], 4), flush=True)
    # fixed-budget design curves from the pooled detection-vs-length measurements
    phis = np.array(PHIS)
    for view, o in out.items():
        pd = [o['thinning'][str(p)]['p_det'] for p in PHIS]
        f = interpolated(phis, pd)
        o['budget'] = {f'order{m}': fixed_budget(m, [2, 3, 4, 6, 8, 12, 16], f) for m in (2, 3)}
    (OUT / 'cd20-design.json').write_text(json.dumps(out, indent=1))


def partitions(k):
    seen = set()
    for assign in itertools.product(range(k), repeat=C.K):
        if len(set(assign)) != k:
            continue
        seen.add(tuple(sorted(tuple(l for l in range(C.K) if assign[l] == g) for g in range(k))))
    return [list(map(list, g)) for g in sorted(seen)]


if __name__ == '__main__':
    main(Path(sys.argv[1]))
