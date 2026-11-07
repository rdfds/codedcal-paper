"""Shape of CD20 assemblies from the angles of connected three-label triplets.

For each view, the angle j-i-k between a molecule's two labelled assembly
partners (within 30 nm, three distinct labels) is profiled after removing every
pair-level contribution. Reference profiles come from simulated straight and flexible
chains, hexameric rings and compact clusters with the same label size (5 nm per
molecule, from the Nup96 fit), labelling efficiency and precision. Shapes are
compared on normalized profiles and on a linearity index: the share of connected
triplets at angles above 150 degrees (straight chains place a molecule between
two partners).

Usage: python shape_cd20.py CACHE_DIR
"""
import json
import sys
from pathlib import Path
import numpy as np
import cd20 as C
from simulate import acquire, complexes
from verify_angles import LINK_SD, oracle_profile

OUT = Path(__file__).resolve().parent / 'results'
B = 400
MODELS = {
    'straight chains': dict(fractions=[0.2, 0.2, 0.15, 0.15, 0.15, 0.15], shape='chain'),
    'hexameric rings': dict(fractions=[0.4, 0.2, 0, 0, 0, 0.4], shape='ring'),
    'compact clusters': dict(fractions=[0.3, 0.2, 0, 0.15, 0.2, 0.15], shape='blob'),
    'flexible chains, 30 deg': dict(fractions=[0.2, 0.2, 0.15, 0.15, 0.15, 0.15], shape='flexible', bend=np.deg2rad(30)),
    'flexible chains, 60 deg': dict(fractions=[0.2, 0.2, 0.15, 0.15, 0.15, 0.15], shape='flexible', bend=np.deg2rad(60)),
    'flexible chains, 90 deg': dict(fractions=[0.2, 0.2, 0.15, 0.15, 0.15, 0.15], shape='flexible', bend=np.deg2rad(90)),
}
centres = (C.ANGLE_EDGES[:-1] + C.ANGLE_EDGES[1:]) / 2
LINEAR = centres > np.deg2rad(150)


def model_profiles(rng, reps=3):
    out = {}
    for name, m in MODELS.items():
        prof = 0
        for _ in range(reps):
            sub, cid = complexes(rng, 10000.0, 150, m['fractions'], 13.5, shape=m['shape'], bend=m.get('bend', 0.0))
            d = acquire(rng, sub, cid, C.K, efficiency=0.4, link_sd=LINK_SD, locs_per_molecule=1, loc_sd=0.0)
            xy = d['site'] + rng.normal(0, 0.5, d['site'].shape)
            prof = prof + oracle_profile(xy, d['site_label'], d['site_complex'])
        out[name] = (prof / prof.sum()).tolist()
    return out


def main(cache):
    rng = np.random.default_rng(8)
    models = model_profiles(rng)
    out = dict(angle_edges=C.ANGLE_EDGES.tolist(), models=models,
               model_linearity={k: float(np.sum(np.array(v)[LINEAR])) for k, v in models.items()}, views={})
    for view in C.VIEWS:
        mols = C.load_resi(cache, view)
        xy = np.concatenate([m['xy'] for m in mols])
        lab = np.concatenate([np.full(len(m['xy']), r) for r, m in enumerate(mols)])
        res = C.point_pair_g(xy, lab, rng)
        est, boot = C.triplet_angles(xy, lab, res, rng, B=B)
        tot, tot_b = est.sum(), boot.sum(1)
        norm = est / tot
        norm_b = boot / tot_b[:, None]
        lin_b = norm_b[:, LINEAR].sum(1)
        se = norm_b.std(0)
        fit = {k: float((((norm - np.array(v)) / np.maximum(se, 1e-6)) ** 2).sum()) for k, v in models.items()}
        out['views'][view] = dict(treated=C.TREATED[view], profile=est.tolist(),
                                  profile_ci=np.percentile(boot, [2.5, 97.5], axis=0).tolist(),
                                  normalized=norm.tolist(), normalized_se=se.tolist(),
                                  linearity=float(norm[LINEAR].sum()),
                                  linearity_ci=np.percentile(lin_b, [2.5, 97.5]).tolist(),
                                  chi2_vs_models=fit, best_model=min(fit, key=fit.get))
        v = out['views'][view]
        print(view, 'total', round(float(tot), 4), '| linearity', round(v['linearity'], 3),
              np.round(v['linearity_ci'], 3).tolist(), '| chi2', {k: round(x, 1) for k, x in fit.items()},
              '| best', v['best_model'], flush=True)
    print('model linearity', {k: round(x, 3) for k, x in out['model_linearity'].items()})
    # finer scan of flexible chains: bend angle x link spacing
    scan = {}
    for bend in (40, 50, 60, 70, 80):
        for spacing in (10.0, 13.5, 17.0):
            prof = 0
            for _ in range(3):
                sub, cid = complexes(rng, 10000.0, 150, MODELS['straight chains']['fractions'], spacing,
                                     shape='flexible', bend=np.deg2rad(bend))
                d = acquire(rng, sub, cid, C.K, efficiency=0.4, link_sd=LINK_SD, locs_per_molecule=1, loc_sd=0.0)
                xy = d['site'] + rng.normal(0, 0.5, d['site'].shape)
                prof = prof + oracle_profile(xy, d['site_label'], d['site_complex'])
            prof = prof / prof.sum()
            scan[f'bend{bend}_spacing{spacing:g}'] = {
                v: float((((np.array(o['normalized']) - prof) / np.maximum(o['normalized_se'], 1e-6)) ** 2).sum())
                for v, o in out['views'].items() if o['treated']}
    out['flexible_scan'] = scan
    for v in [k for k, o in out['views'].items() if o['treated']]:
        best = min(scan, key=lambda k: scan[k][v])
        print(v, 'best flexible chain', best, round(scan[best][v], 1))
    (OUT / 'cd20-shape.json').write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main(Path(sys.argv[1]))
