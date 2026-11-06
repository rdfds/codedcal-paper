"""Label permanence: is any molecule recorded in more than one round?

For every published RESI molecule, the nearest molecule in each other round is
found. A molecule seen in two rounds would give a cross-round neighbour at the
combined precision and registration error (~1-2 nm). The observed fraction of
cross-round nearest neighbours within CLOSE_NM is compared with the same
fraction after displacing the other round by 1.5-3 um (no shared molecules),
with and without removing the per-round translations estimated from the raw
localizations (results/cd20-derived.json).

Usage: python permanence_cd20.py CACHE_DIR
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import cd20 as C

OUT = Path(__file__).resolve().parent / 'results'
CLOSE_NM = 2.0
EDGES = np.arange(0, 20.5, 1.0)


def nn_hist(mols, shift=None):
    h = np.zeros(len(EDGES) - 1); close = 0; total = 0
    for a, A in enumerate(mols):
        for b, Bm in enumerate(mols):
            if a == b:
                continue
            pts = Bm if shift is None else Bm + shift[b] - shift[a]
            d, _ = cKDTree(pts).query(A)
            h += np.histogram(d, EDGES)[0]
            close += (d < CLOSE_NM).sum(); total += len(d)
    return h, close / total


def main(cache):
    rng = np.random.default_rng(5)
    derived = json.loads((OUT / 'cd20-derived.json').read_text())
    out = {}
    for view in C.VIEWS:
        raw = [m['xy'] for m in C.load_resi(cache, view)]
        t = np.array(derived[view]['round_translation_nm'])
        corrected = [xy - t[r] for r, xy in enumerate(raw)]
        res = {}
        for name, mols in (('as_published', raw), ('translation_corrected', corrected)):
            h, frac = nn_hist(mols)
            null = [nn_hist(mols, s) for s in C.null_shifts(rng, 4)]
            hn = np.mean([x[0] for x in null], 0); fn = float(np.mean([x[1] for x in null]))
            # ratio only where the null expects at least one neighbour in the bin
            ratio = [round(float(a / b), 3) if b >= 1 else None for a, b in zip(h, hn)]
            res[name] = dict(close_fraction=float(frac), null_close_fraction=fn,
                             excess_close_fraction=float(frac - fn),
                             observed_by_nm=h.astype(int).tolist(), null_by_nm=np.round(hn, 2).tolist(),
                             ratio_by_nm=ratio)
        out[view] = res
        print(view, {k: (round(v['close_fraction'], 4), round(v['null_close_fraction'], 4))
                     for k, v in res.items()}, flush=True)
    (OUT / 'cd20-permanence.json').write_text(json.dumps(dict(close_nm=CLOSE_NM, edges=EDGES.tolist(),
                                                              views=out), indent=1))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
