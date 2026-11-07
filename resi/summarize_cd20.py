"""Summaries from saved per-view pair counts (no recounting).

- excess labelled neighbours within 30 nm (pooled, per label pair, merged labels)
- per-round translations (from results/cd20-registration.json, which checks
  each view by injecting a known shift)

Usage: python summarize_cd20.py CACHE_DIR
"""
import itertools
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import cd20 as C

OUT = Path(__file__).resolve().parent / 'results'
B = 400


def local_density(cache, view, annulus=(300.0, 600.0)):
    """Pair-weighted density of RESI molecules (all rounds) on the plateau annulus,
    per nm^2: the density scale implied by plateau normalization."""
    xy = np.concatenate([r['xy'] for r in C.load_resi(cache, view)])
    t = cKDTree(xy)
    cnt = t.count_neighbors(t, np.array(annulus))
    return float(np.diff(cnt)[0] / len(xy) / (np.pi * (annulus[1] ** 2 - annulus[0] ** 2)))


def merged(N, n, groups):
    """Pair counts after merging labels into groups (list of label lists)."""
    k = len(groups)
    Nm = np.zeros(N.shape[:1] + (k, k) + N.shape[3:])
    nm = np.zeros((n.shape[0], k))
    for i, gi in enumerate(groups):
        nm[:, i] = n[:, gi].sum(1)
        for j, gj in enumerate(groups):
            Nm[:, i, j] = N[:, gi][:, :, gj].sum((1, 2))
    return Nm, nm


def pooled_cross(N, n):
    k = N.shape[1]
    Ns, ns = N.sum(0), n.sum(0)
    off = ~np.eye(k, dtype=bool)
    return Ns[off].sum(0) / np.outer(ns, ns)[off].sum()


def excess_stat(N, n, rho, null_density, plateau_ratio):
    g = C.g_of_r(pooled_cross(N, n), null_density, plateau_ratio)
    return C.excess_neighbours(g, rho)


def boot(fn, N, n, rng):
    nb = N.shape[0]
    return np.array([fn(N[k], n[k]) for k in (rng.integers(0, nb, nb) for _ in range(B))])


def ci(x):
    lo, hi = np.percentile(x, [2.5, 97.5])
    return [round(float(lo), 4), round(float(hi), 4)]


def main(cache):
    rng = np.random.default_rng(1)
    registration = json.loads((OUT / 'cd20-registration.json').read_text())
    out = {}
    for view in C.VIEWS:
        d = np.load(OUT / f'cd20-{view}.npz')
        N, n = d['N'], d['n']
        nd, pr = float(d['null_density']), float(d['plateau_ratio'])
        rho = local_density(cache, view)
        stat = lambda a, b: excess_stat(a, b, rho, nd, pr)
        pooled = stat(N, n)
        pooled_b = boot(stat, N, n, rng)
        pairs = {}
        for a, b in itertools.combinations(range(C.K), 2):
            Np, npp = merged(N, n, [[a], [b]])
            pairs[f'R{a+1}R{b+1}'] = [round(stat(Np, npp), 4), ci(boot(stat, Np, npp, rng))]
        merges = {}
        for k in (3, 2):
            vals = []
            for groups in set_partitions(k):
                Nm, nm = merged(N, n, groups)
                vals.append(stat(Nm, nm))
            merges[f'K={k}'] = [round(float(np.mean(vals)), 4), round(float(np.min(vals)), 4),
                                round(float(np.max(vals)), 4), len(vals)]
        reg = registration[view]
        out[view] = dict(treated=C.TREATED[view], local_density_per_um2=round(rho * 1e6, 1),
                         excess_neighbours_30nm=[round(pooled, 4), ci(pooled_b)],
                         excess_by_label_pair=pairs, excess_merged_labels=merges,
                         round_translation_nm=reg['round_translation_nm'],
                         max_round_shift_nm=round(reg['max_round_shift_nm'], 2),
                         translation_fit_rms_nm=round(reg['translation_fit_rms_nm'], 3),
                         registration_reliable=reg['reliable'])
        print(view, json.dumps(out[view]), flush=True)
    (OUT / 'cd20-derived.json').write_text(json.dumps(out, indent=1))


def set_partitions(k):
    """Partitions of the 4 labels into k nonempty groups."""
    labels = list(range(C.K))
    seen = set()
    for assign in itertools.product(range(k), repeat=C.K):
        if len(set(assign)) != k:
            continue
        groups = tuple(sorted(tuple(l for l in labels if assign[l] == g) for g in range(k)))
        seen.add(groups)
    return [list(map(list, g)) for g in sorted(seen)]


if __name__ == '__main__':
    main(Path(sys.argv[1]))
