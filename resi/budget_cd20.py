"""Fixed-budget design curves for several minimum-localization rules.

Detection probability versus relative round length phi is measured by thinning
the raw localizations of each published RESI molecule (see design_cd20.py) for
each threshold; the relative pair (order 2) and triplet (order 3) information
of K equal rounds in the total time of four full rounds follows from design.py.
Detection is relative to the published molecules, so lower thresholds do not
include molecules the original pipeline missed, nor any false detections.

Usage: python budget_cd20.py CACHE_DIR
"""
import json
import sys
from pathlib import Path
import numpy as np
import cd20 as C
from design import fixed_budget, interpolated
from design_cd20 import detect, molecules_from_rounds

OUT = Path(__file__).resolve().parent / 'results'
PHIS = [1.0, 0.5, 0.25, 0.125, 0.0625]
THRESHOLDS = [16, 8, 4]
KS = [2, 3, 4, 6, 8, 12, 16]


def main(cache):
    out = {}
    for view in C.VIEWS:
        rounds = molecules_from_rounds(cache, view)
        n_resi = sum(d['n_mol'] for d in rounds)
        row = {}
        for t in THRESHOLDS:
            pdet = [len(detect(rounds, phi, t)[0]) / n_resi for phi in PHIS]
            f = interpolated(PHIS, pdet)
            row[f'min{t}'] = dict(p_det=dict(zip(map(str, PHIS), np.round(pdet, 4).tolist())),
                                  order2=fixed_budget(2, KS, f), order3=fixed_budget(3, KS, f))
            best2 = max(row[f'min{t}']['order2'], key=row[f'min{t}']['order2'].get)
            best3 = max(row[f'min{t}']['order3'], key=row[f'min{t}']['order3'].get)
            row[f'min{t}'].update(best_K_order2=best2, best_K_order3=best3)
            print(view, f'min{t}', 'p_det', np.round(pdet, 3).tolist(), 'best K pairs', best2,
                  'triplets', best3, 'gain3', round(row[f'min{t}']['order3'][best3], 2), flush=True)
        out[view] = row
    (OUT / 'cd20-budget.json').write_text(json.dumps(dict(phis=PHIS, K=KS, views=out), indent=1))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
