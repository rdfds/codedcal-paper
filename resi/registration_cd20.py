"""Fiducial-free round registration from cross-label pairs, with injection checks.

For each view (raw round localizations), offsets between rounds are estimated
with crosslabel.registration_offsets (plane-corrected centroid) and solved for
per-round translations. Each view is then re-analysed after adding a known
shift to one round; the change in that round's translation must recover the
injected shift. A view's offsets are reported as reliable only if recovery is
within RECOVERY_TOLERANCE in both axes.

Usage: python registration_cd20.py CACHE_DIR
"""
import json
import sys
from pathlib import Path
import numpy as np
import cd20 as C

OUT = Path(__file__).resolve().parent / 'results'
INJECT = np.array([2.0, -1.5])     # nm, added to round 3
INJECT_ROUND = 2
RECOVERY_TOLERANCE = 0.25


def load_npc(cache):
    xy, lab = [], []
    for r in range(1, C.K + 1):
        d = np.load(Path(cache) / f'npc-locs-R{r}.npz')
        keep = d['lp'] <= C.LP_MAX
        xy.append(np.column_stack([d['x'][keep], d['y'][keep]]).astype(float)); lab.append(np.full(keep.sum(), r - 1))
    return np.concatenate(xy), np.concatenate(lab)


def translations(xy, lab):
    off, se, cnt = C.registration_offsets(xy, lab)
    t, rms = C.label_translations(off)
    return off, se, t, rms


def main(cache):
    out = {}
    views = C.VIEWS + ['npc']
    for view in views:
        if view == 'npc':
            xy, lab = load_npc(cache)
        else:
            d = C.load(cache, view); xy, lab = d['xy'], d['label']
        off, se, t, rms = translations(xy, lab)
        moved = xy.copy(); moved[lab == INJECT_ROUND] += INJECT
        _, _, t2, _ = translations(moved, lab)
        others = [r for r in range(C.K) if r != INJECT_ROUND]
        rel = lambda tt: tt[INJECT_ROUND] - tt[others].mean(0)
        recovered = rel(t2) - rel(t)
        ratio = recovered / INJECT
        reliable = bool(np.all(np.abs(ratio - 1) <= RECOVERY_TOLERANCE))
        out[view] = dict(treated=C.TREATED.get(view), round_translation_nm=np.round(t, 3).tolist(),
                         max_round_shift_nm=float(np.max(np.hypot(*t.T))), translation_fit_rms_nm=rms,
                         pair_offsets_nm={f'R{a+1}-R{b+1}': [np.round(off[a, b], 3).tolist(), np.round(se[a, b], 3).tolist()]
                                          for a in range(C.K) for b in range(a + 1, C.K)},
                         injected_nm=INJECT.tolist(), recovered_nm=np.round(recovered, 3).tolist(),
                         recovery_ratio=np.round(ratio, 3).tolist(), reliable=reliable)
        print(view, '| max shift', round(out[view]['max_round_shift_nm'], 2), 'nm | additivity rms', round(rms, 2),
              '| injection recovered', np.round(recovered, 2).tolist(), '| reliable', reliable, flush=True)
    (OUT / 'cd20-registration.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
