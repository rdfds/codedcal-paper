"""Independent Nup96 NPC validation on data from another laboratory and microscope.

Zaza et al., Nat Commun 16:4991 (2025), Fig. 3 RESI: U2OS Nup96-mEGFP labelled
with a pool of anti-GFP nanobodies carrying four docking sequences, imaged
round by round on a spinning-disk confocal with optical photon reassignment
(figshare 10.6084/m9.figshare.28741373). The published RESI molecule table
carries each molecule's round (resi_channel_id). Without curated picks, pores
are found as connected components of molecules (single linkage, LINK_NM);
components larger than an NPC (merged neighbours) are set aside. With about a
third of copies observed, shorter links split rings at empty corners; pore-level
estimates stabilize from ~80 nm (scan recorded in the output). The analysis then follows
npc_validation.py: label linkage spread from the corner structure, geometry
factors, and labelling-and-detection efficiency by every route.

Usage: python zaza_npc.py CACHE_DIR
"""
import json
import sys
from pathlib import Path
import h5py
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
import npc as P
from npc_validation import NN_EDGES, RADIUS, estimate, nn_hist

OUT = Path(__file__).resolve().parent / 'results'
URL = 'https://ndownloader.figshare.com/files/53497595'
MEMBER = 'Data_submission/Fig 3/RESI/NUP96_resi.hdf5'
PIXEL_NM = 108.0         # 60x objective, 4x SoRa magnification, 4x4 binning of 6.5-um pixels
LINK_NM = 80.0           # links copies of one pore across empty corners (scan in output)
MAX_EXTENT_NM = 170.0    # an NPC spans ~100 nm plus label spread


def fetch(cache):
    dest = Path(cache) / 'zaza-NUP96_resi.hdf5'
    if not dest.exists():
        from remote_zip import open_remote_zip
        dest.write_bytes(open_remote_zip(URL).read(MEMBER))
    return dest


def load(cache):
    with h5py.File(fetch(cache)) as f:
        a = f['locs'][:]
    xy = np.column_stack([a['x'], a['y']]).astype(float) * PIXEL_NM
    lab = a['resi_channel_id'].astype(int)
    pairs = cKDTree(xy).query_pairs(LINK_NM, output_type='ndarray')
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(xy), len(xy)))
    _, comp = connected_components(graph, directed=False)
    keep = np.zeros(len(xy), bool)
    sizes = np.bincount(comp)
    for c in np.flatnonzero(sizes >= 2):
        m = comp == c
        p = xy[m]
        if np.ptp(p, axis=0).max() <= MAX_EXTENT_NM:
            keep[m] = True
    ring = [np.hypot(*(xy[comp == c] - xy[comp == c].mean(0)).T).mean()
            for c in np.unique(comp[keep]) if (comp == c).sum() >= 8]
    info = dict(molecules=int(len(xy)), assigned=int(keep.sum()), pores=int(len(np.unique(comp[keep]))),
                single_molecule_components=int((sizes == 1).sum()),
                mean_distance_to_pore_centre_nm=float(np.mean(ring)))
    return xy[keep], lab[keep], comp[keep], info


def link_scan(cache):
    global LINK_NM
    keep, out = LINK_NM, {}
    for link in (45.0, 60.0, 70.0, 80.0, 90.0, 100.0):
        LINK_NM = link
        xy, lab, pid, info = load(cache)
        pi = np.bincount(lab, minlength=P.K) / len(lab)
        st = P.stoichiometry(P.pore_counts(lab, pid), pi)
        out[f'{link:g}'] = dict(pores=info['pores'], assigned=info['assigned'],
                                **{k: float(v) for k, v in st.items()})
    LINK_NM = keep
    return out


def main(cache):
    rng = np.random.default_rng(9)
    scan = link_scan(cache)
    xy, lab, pid, info = load(cache)
    pi = np.bincount(lab, minlength=P.K) / len(lab)
    obs_h = nn_hist(xy, lab, pid)
    fits = {}
    for sd in (2.0, 3.0, 4.0, 5.0, 6.0, 7.0):
        sx, sl, sp = P.simulate(rng, 600, 0.4, pi, link_sd=sd)
        sim_h = nn_hist(sx, sl, sp)
        fits[sd] = float((sim_h - obs_h) @ (sim_h - obs_h))
    link_sd = min(fits, key=fits.get)
    a1, a2 = P.geometry_factors(rng, RADIUS, link_sd=link_sd)
    data = estimate(xy, lab, pid, pi, a1, a2, rng)
    print('info', info)
    print('label fractions', pi.round(3).tolist(), '| link sd', link_sd, '| a1, a2', round(a1, 3), round(a2, 3))
    print('data', {k: round(v, 3) for k, v in data.items() if k != 'ci'})
    print('data CIs', {k: np.round(v, 3).tolist() for k, v in data['ci'].items()})
    out = dict(source='Zaza et al., Nat Commun 16:4991 (2025); figshare 10.6084/m9.figshare.28741373',
               pixel_nm=PIXEL_NM, link_nm=LINK_NM, max_extent_nm=MAX_EXTENT_NM, info=info,
               label_fractions=pi.round(4).tolist(), radius_nm=RADIUS, link_sd_fit=link_sd,
               link_sd_misfit=fits, a1=a1, a2=a2, observed_nn_hist=obs_h.tolist(),
               nn_edges=NN_EDGES.tolist(), data=data, link_scan=scan)
    (OUT / 'npc-zaza.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
