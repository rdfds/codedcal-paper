"""Recover a known dimer fraction with the full CD20 pipeline on simulated fields.

Simulated fields use the CD20 pipeline unchanged: block pair counts, shifted
null, plateau normalization, block bootstrap and the monomer/dimer fit. Rounds
have unequal localization yields (label-dependent counts), as in the real data.
"""
import json
import sys
from pathlib import Path
import numpy as np
import cd20 as C
from model import bootstrap_fit, fit_dimer
from simulate import acquire, complexes

SIDE = 10000.0
SUBUNIT_DENSITY = 100.0        # per um^2
EFFICIENCY = 0.5
LOC_SD = 4.0
LINK_SD = 4.0
GAIN = [1.0, 0.7, 1.1, 1.5]


def one(f_dimer_subunits, seed):
    rng = np.random.default_rng(seed)
    # fraction of complexes that are dimers giving the requested subunit fraction
    q = f_dimer_subunits / (2 - f_dimer_subunits)
    sub, cid = complexes(rng, SIDE, SUBUNIT_DENSITY, [1 - q, q], 13.5)
    data = acquire(rng, sub, cid, C.K, efficiency=EFFICIENCY, link_sd=LINK_SD,
                   locs_per_molecule=30, round_gain=GAIN, loc_sd=LOC_SD)
    xy, lab = data['xy'], data['label']
    block = C.blocks(xy)
    N, n = C.block_pair_counts(xy, lab, block=block)
    null_density, plateau_ratio = C.null_normalization(xy, lab, rng, n_shift=4)
    _, pooled = C.cross_g(N, n)
    g = C.g_of_r(pooled, null_density, plateau_ratio)
    boot = C.bootstrap(lambda a, b: C.g_of_r(C.cross_g(a, b)[1], null_density, plateau_ratio),
                       N, n, rng, B=100)
    se = boot.std(0)
    fit = fit_dimer(C.EDGES, g, se, 2 * LOC_SD ** 2)
    bf = bootstrap_fit(C.EDGES, boot, se, 2 * LOC_SD ** 2)
    true_rho = len(sub) / SIDE ** 2
    p = fit['params']
    f_lo, f_hi = np.percentile(bf[:, 0] * true_rho, [2.5, 97.5])
    d_lo, d_hi = np.percentile(bf[:, 1], [2.5, 97.5])
    realized = float(np.mean(np.bincount(cid)[cid] == 2))
    xs = C.excess_neighbours(g, true_rho)
    xs_lo, xs_hi = np.percentile(C.excess_neighbours(boot, true_rho), [2.5, 97.5])
    return dict(true_f=f_dimer_subunits, realized_f=realized,
                f_hat=p['amplitude_nm2'] * true_rho, f_ci=[f_lo, f_hi],
                f_covered=bool(f_lo <= realized <= f_hi),
                excess=xs, excess_ci=[xs_lo, xs_hi], excess_covered=bool(xs_lo <= realized <= xs_hi),
                d_hat=p['d_nm'], d_ci=[d_lo, d_hi], link_hat=p['sigma_link_nm'],
                baseline=p.get('baseline'))


if __name__ == '__main__':
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / 'results'
    out.mkdir(exist_ok=True)
    reps = {0.0: 8, 0.3: 4, 0.47: 4, 0.64: 4}
    rows = [one(f, 50 + 10 * k + rep) for k, (f, n) in enumerate(reps.items()) for rep in range(n)]
    (out / 'verify_model.json').write_text(json.dumps(rows, indent=1))
    for r in rows:
        print({k: (round(v, 3) if isinstance(v, float) else
                   [round(x, 3) for x in v] if isinstance(v, list) else v) for k, v in r.items()})
    print('fit coverage:', np.mean([r['f_covered'] for r in rows]),
          '| excess-neighbour coverage:', np.mean([r['excess_covered'] for r in rows]))
