"""Experimental bootstrap-selection sensitivity; not an error probability."""
import numpy as np
from .output_risk import canonical_maps,catalog_distance

def selection_weights(profile,repetitions=500,seed=61021,tie_tolerance=1e-6):
    residual=np.asarray(profile['residual']);energy=np.asarray(profile['energy']);counts=np.asarray(profile['block_counts']);use=(counts>0)&(energy.max(axis=0)>0);n=int(use.sum())
    if n<12 or profile['coverage']<.8:return dict(abstained=True,reason='Insufficient common supported blocks')
    residual=residual[:,use];energy=energy[:,use];rng=np.random.default_rng(seed);index=rng.integers(0,n,(repetitions,n));draws=np.stack([residual[:,i].sum(axis=1)/np.maximum(energy[:,i].sum(axis=1),1e-20) for i in index],axis=1);ties=draws<=draws.min(axis=0)+tie_tolerance;weights=(ties/ties.sum(axis=0)).mean(axis=1);assert abs(weights.sum()-1)<1e-12
    return dict(abstained=False,weights=weights.tolist(),tie_tolerance=tie_tolerance,multiple_tie_draws=int((ties.sum(axis=0)>1).sum()),seed=seed,repetitions=repetitions,posterior_probability=False)

def weighted_output_risk(catalogs,selection,incumbent,maps):
    if selection['abstained']:return dict(abstained=True,reason=selection['reason'])
    base=catalogs[incumbent];coding=np.sum(~np.isin(base['words'],base.get('negative_ids',[])))
    if coding<10:return dict(abstained=True,reason='Insufficient incumbent coding calls')
    weights=np.asarray(selection['weights']);changes=[catalog_distance(base,c) for c in catalogs];identity=float(weights @ np.array([c['identity_loss'] for c in changes]));tv=float(weights @ np.array([c['composition_tv'] for c in changes]));relative=canonical_maps(maps);scale=np.array([12.,48.,48.]);delta=(relative-relative[incumbent])/scale;geometry=float(np.sqrt(weights @ np.square(delta).sum(axis=2).mean(axis=1)))
    return dict(abstained=False,identity_risk=identity,composition_risk=tv,flagged=identity>=.05 or tv>.05,weighted_geometric_RMS=geometry,selection_weights=weights.tolist(),molecular_correctness_certified=False,posterior_RNA_error_probability=False)
