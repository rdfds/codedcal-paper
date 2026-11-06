"""Verify the exact label-draw identities on fixed simulated configurations.

For each configuration, molecule labels are redrawn many times while positions
and counts stay fixed. The Monte Carlo mean of each statistic is compared with
its exact conditional expectation, and the Monte Carlo variance of the cross
statistic with its exact variance.
"""
import json
import sys
from pathlib import Path
import numpy as np
from labels import (decompose, exact_cross_variance, label_moments, oracle_pairs,
                    oracle_star, pair_counts, split_pairs, star_triplets)
from simulate import acquire, complexes, relabel

EDGES = np.arange(0, 41, 5.0)
RADIUS = 20.0
DRAWS = 400

CASES = [
    dict(name='uniform K=4, monomer/dimer', K=4, pi=None, fractions=[0.7, 0.3]),
    dict(name='uniform K=2, with trimers', K=2, pi=None, fractions=[0.5, 0.3, 0.2]),
    dict(name='non-uniform K=4, with trimers', K=4, pi=[0.4, 0.3, 0.2, 0.1], fractions=[0.5, 0.3, 0.2]),
    dict(name='uniform K=8, with trimers', K=8, pi=None, fractions=[0.5, 0.3, 0.2]),
]


def run(case, seed):
    rng = np.random.default_rng(seed)
    K = case['K']
    pi = np.full(K, 1 / K) if case['pi'] is None else np.asarray(case['pi'])
    sub, cid = complexes(rng, side=1500, density=120, fractions=case['fractions'], spacing=13.5)
    data = acquire(rng, sub, cid, K, efficiency=0.6, pi=pi, link_sd=3.0,
                   locs_per_molecule=12, loc_sd=4.0)
    xy, mol = data['xy'], data['mol']
    S, W, pairs = oracle_pairs(xy, mol, EDGES)
    W3 = oracle_star(xy, mol, RADIUS)
    m = label_moments(pi)
    same, cross, star, Shat, What = [], [], [], [], []
    for _ in range(DRAWS):
        lab = relabel(rng, data, K, pi)
        N = pair_counts(xy, lab, EDGES, K)
        sm, cr = split_pairs(N)
        dec = decompose(N, pi)
        same.append(sm); cross.append(cr); Shat.append(dec['S']); What.append(dec['W'])
        star.append(star_triplets(xy, lab, RADIUS, K))
    same, cross, star = map(np.array, (same, cross, star))
    Shat, What = np.array(Shat), np.array(What)
    var_exact = exact_cross_variance(pairs, pi, len(EDGES) - 1, int(mol.max()) + 1)

    def zscore(mc, exact):
        se = mc.std(0, ddof=1) / np.sqrt(len(mc))
        return (mc.mean(0) - exact) / np.where(se > 0, se, np.inf)

    star_exact = m['distinct3'] * W3
    if m['distinct3'] < 1e-12:  # K = 2: three distinct labels are impossible
        assert not star.any()
        z_star, star_ratio = 0.0, None
    else:
        z_star = float(zscore(star[:, None], np.array([star_exact]))[0])
        star_ratio = float(star.mean() / star_exact)
    return dict(
        case=case['name'], localizations=len(xy), molecules=int(mol.max()) + 1,
        z_same=zscore(same, S + m['s2'] * W).tolist(),
        z_cross=zscore(cross, m['cross'] * W).tolist(),
        z_S=zscore(Shat, S).tolist(), z_W=zscore(What, W).tolist(),
        z_star=z_star,
        var_ratio=(cross.var(0, ddof=1) / np.where(var_exact > 0, var_exact, np.nan)).tolist(),
        star_mean_over_exact=star_ratio,
    )


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / 'results'
    out.mkdir(exist_ok=True)
    results = [run(c, 1000 + k) for k, c in enumerate(CASES)]
    zs = np.concatenate([np.r_[r['z_same'], r['z_cross'], r['z_S'], r['z_W'], r['z_star']]
                         for r in results])
    zs = zs[np.isfinite(zs)]
    vr = np.concatenate([r['var_ratio'] for r in results])
    vr = vr[np.isfinite(vr)]
    summary = dict(draws=DRAWS, n_z=len(zs), max_abs_z=float(np.abs(zs).max()),
                   frac_abs_z_gt_3=float((np.abs(zs) > 3).mean()),
                   var_ratio_min=float(vr.min()), var_ratio_max=float(vr.max()),
                   var_ratio_median=float(np.median(vr)))
    (out / 'verify_exact.json').write_text(json.dumps(dict(summary=summary, cases=results), indent=2))
    print(json.dumps(summary, indent=2))
    for r in results:
        ratio = r['star_mean_over_exact']
        print(r['case'], '| locs', r['localizations'],
              '| star mean/exact', 'n/a (K=2)' if ratio is None else round(ratio, 4),
              '| var ratio', np.round(r['var_ratio'], 3).tolist())
