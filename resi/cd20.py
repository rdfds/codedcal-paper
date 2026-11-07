"""Loading of the public CD20 RESI acquisitions and the cross-label statistics
used on them, with the four labels (rounds R1-R4) fixed.

Generic statistics live in the crosslabel package; functions whose behaviour
depends on the number of labels are bound to K = 4 here, so that merged-label
analyses keep four-direction null shifts and 4 x 4 label-pair arrays.
"""
from functools import partial
from pathlib import Path
import numpy as np
from crosslabel import correlation as _c
from crosslabel.correlation import (ANGLE_EDGES, BLOCK, EDGES, NULL_SMOOTH, PLATEAU, PLATEAU_FRACTION,
                                    TRIPLET_ANNULUS, TRIPLET_RADIUS, baseline_excess, blocks, bootstrap,
                                    cross_g, decomposition, excess_neighbours, excess_with_ci, g_of_r,
                                    label_translations, plateau_density, triplet_angles, triplet_m2)

VIEWS = ['fig4untreated', 'edf10rep1', 'edf10rep2', 'rtx_rep1', 'rtx_rep2', 'rtx_rep3']
TREATED = {'fig4untreated': False, 'edf10rep1': False, 'edf10rep2': False,
           'rtx_rep1': True, 'rtx_rep2': True, 'rtx_rep3': True}
K = 4
LP_MAX = 8.0                      # nm, per-localization precision gate

block_pair_counts = partial(_c.block_pair_counts, n_labels=K)
cross_annuli = partial(_c.cross_annuli, n_labels=K)
null_shifts = partial(_c.null_shifts, n_labels=K)
null_normalization = partial(_c.null_normalization, n_labels=K)
point_pair_g = partial(_c.point_pair_g, n_labels=K)
registration_offsets = partial(_c.registration_offsets, n_labels=K)


def load(cache, view, frame_fraction=1.0, lp_max=LP_MAX):
    """Concatenated localizations of the four rounds with their round label."""
    xy, lab, frame, lp = [], [], [], []
    for r in range(1, K + 1):
        d = np.load(Path(cache) / f'{view}-locs-R{r}.npz')
        keep = d['lp'] <= lp_max
        if frame_fraction < 1:
            keep &= d['frame'] < frame_fraction * (d['frame'].max() + 1)
        xy.append(np.column_stack([d['x'][keep], d['y'][keep]]).astype(float))
        lab.append(np.full(keep.sum(), r - 1))
        frame.append(d['frame'][keep]); lp.append(d['lp'][keep])
    return dict(xy=np.concatenate(xy), label=np.concatenate(lab),
                frame=np.concatenate(frame), lp=np.concatenate(lp).astype(float))


def load_resi(cache, view):
    out = []
    for r in range(1, K + 1):
        d = np.load(Path(cache) / f'{view}-resi-R{r}.npz')
        out.append(dict(xy=np.column_stack([d['x'], d['y']]), lp=d['lp']))
    return out
