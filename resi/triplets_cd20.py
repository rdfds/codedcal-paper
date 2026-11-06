"""Connected triplet statistic on the published RESI molecule positions.

m2 is the mean number of unordered pairs of labelled partners within 30 nm per
labelled molecule, from three-distinct-label triplets after removing everything
pair correlations can produce. It is zero for any monomer/dimer mixture.

Usage: python triplets_cd20.py CACHE_DIR
"""
import json
import sys
from pathlib import Path
import numpy as np
import cd20 as C

OUT = Path(__file__).resolve().parent / 'results'
B = 400


def main(cache):
    rng = np.random.default_rng(3)
    out = {}
    for view in C.VIEWS:
        mols = C.load_resi(cache, view)
        xy = np.concatenate([m['xy'] for m in mols])
        lab = np.concatenate([np.full(len(m['xy']), r) for r, m in enumerate(mols)])
        res = C.point_pair_g(xy, lab, rng)
        m2, ci, p3 = C.triplet_m2(xy, lab, res, rng, B=B)
        pi = np.bincount(lab, minlength=C.K) / len(lab)
        out[view] = dict(treated=C.TREATED[view], molecules=int(len(xy)),
                         label_fractions=np.round(pi, 3).tolist(), p3=round(p3, 4),
                         m2=round(m2, 4), m2_ci=np.round(ci, 4).tolist())
        print(view, json.dumps(out[view]), flush=True)
    (OUT / 'cd20-triplets.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
