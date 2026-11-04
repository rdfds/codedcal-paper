"""Experimental finite-profile output instability, not a molecular certificate.

Image/codebook fit compatibility uses spatial block resampling on a fixed,
common supported domain. The statistic describes the sampled transform bank;
it neither identifies true RNA labels nor certifies unsearched distortions.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


def canonical_maps(maps: np.ndarray) -> np.ndarray:
    maps = np.asarray(maps)
    if maps.ndim != 3 or maps.shape[-1] != 3 or not np.isfinite(maps).all():
        raise ValueError('Expected finite candidate × frame × ZYX maps')
    if not np.array_equal(maps, np.rint(maps)):
        raise ValueError('This experimental profile requires integer maps')
    maps = maps.astype(np.int64)
    return maps - maps[:, :1]


def profile_losses(images, support, book, maps, region, block_shape):
    """Return residual/energy sums per candidate and fixed spatial block.

    Pull from unchanged inputs; never stitch bootstrap blocks into fake images.
    A spatial site is retained only when *all* searched maps and all frames
    support it. Support is not silently redefined for each correction.
    """
    images = np.asarray(images); support = np.asarray(support, bool)
    book = np.asarray(book, float); maps = canonical_maps(maps)
    if images.ndim != 4 or support.shape != images.shape or book.ndim != 2 or book.shape[0] != len(images) or maps.shape[1] != len(images):
        raise ValueError('Incompatible input, support or book shapes')
    if not np.isfinite(images).all() or np.any(images < 0) or not np.isfinite(book).all() or np.any(book < 0):
        raise ValueError('Expected finite nonnegative conditioned data and book')
    norm = np.sqrt(np.square(book).sum(axis=0))
    if np.any(norm == 0): raise ValueError('Empty codeword')
    normalized = (book / norm).astype(np.float32)
    lo = np.array([a for a,b in region]); hi = np.array([b for a,b in region]); step = np.array(block_shape)
    extent = hi-lo
    if np.any(step <= 0) or np.any(extent <= 0) or np.any(extent % step): raise ValueError('Region must divide into fixed blocks')
    grid = np.stack(np.meshgrid(*[np.arange(a,b) for a,b in region], indexing='ij'),-1).reshape(-1,3)
    blocks = ((grid-lo)//step); dims = extent//step
    labels = np.ravel_multi_index(blocks.T, dims); count = int(np.prod(dims)); common = np.ones(len(grid),bool)
    for m in maps:
        for f, shift in enumerate(m):
            src = grid-shift
            inside = np.all((src >= 0)&(src < np.array(images.shape[1:])),axis=1)
            safe = src.clip(0,np.array(images.shape[1:])-1)
            common &= inside & support[f,safe[:,0],safe[:,1],safe[:,2]]
    coverage = float(common.mean()); grid=grid[common]; labels=labels[common]
    pixel_counts=np.bincount(labels,minlength=count); residual=[];energy=[]
    for m in maps:
        values=np.stack([images[f,(grid-shift)[:,0],(grid-shift)[:,1],(grid-shift)[:,2]] for f,shift in enumerate(m)]).astype(np.float32)
        e=np.square(values.astype(float)).sum(axis=0)
        # Code chunks bound memory even for a complete 998-word culture book.
        best=np.zeros(len(grid),np.float32)
        for start in range(0,book.shape[1],128):
            best=np.maximum(best,(normalized[:,start:start+128].T @ values).max(axis=0))
        r=np.maximum(e-np.square(best.astype(float)),0)
        residual.append(np.bincount(labels,weights=r,minlength=count));energy.append(np.bincount(labels,weights=e,minlength=count))
    return dict(residual=np.asarray(residual),energy=np.asarray(energy),block_counts=pixel_counts,coverage=coverage,canonical_maps=maps)


def compatibility(profile, repetitions=500, seed=61021, tolerance=.002, min_blocks=12, min_coverage=.8):
    residual=np.asarray(profile['residual']);energy=np.asarray(profile['energy']);counts=np.asarray(profile['block_counts']);size=len(counts)
    usable=(counts>0)&(energy.max(axis=0)>0);n=int(usable.sum())
    if n<min_blocks or profile['coverage']<min_coverage:
        return dict(abstained=True,reason='Insufficient common supported spatial blocks',usable_blocks=n,coverage=profile['coverage'])
    residual=residual[:,usable];energy=energy[:,usable]
    score=residual.sum(axis=1)/np.maximum(energy.sum(axis=1),1e-20);best=int(np.argmin(score))
    rng=np.random.default_rng(seed);indices=rng.integers(0,n,(repetitions,n))
    weights=np.stack([np.bincount(i,minlength=n) for i in indices]);boot=(residual @ weights.T)/np.maximum(energy @ weights.T,1e-20)
    difference=boot-boot[best];lower=np.quantile(difference,.025,axis=1);upper=np.quantile(difference,.975,axis=1)
    compatible=lower<=tolerance;winner=boot.argmin(axis=0);freq=np.bincount(winner,minlength=len(score))/repetitions
    return dict(abstained=False,coverage=profile['coverage'],usable_blocks=n,losses=score.tolist(),best=best,compatible=compatible.tolist(),paired_difference_lower=lower.tolist(),paired_difference_upper=upper.tolist(),bootstrap_winner_frequency=freq.tolist(),tolerance=tolerance,seed=seed,bootstrap_repetitions=repetitions,coverage_guarantee=False)


def catalog_distance(a,b,radius=3.):
    """Identity and composition disagreements; no truth labels supplied."""
    ax=np.asarray(a['xyz'],float).reshape(-1,3);bx=np.asarray(b['xyz'],float).reshape(-1,3);aw=np.asarray(a['words'],int);bw=np.asarray(b['words'],int)
    if len(ax)!=len(aw) or len(bx)!=len(bw): raise ValueError('Mismatched catalog')
    if not len(aw) and not len(bw):return dict(identity_loss=0.,composition_tv=0.,empty_both=True)
    if not len(aw) or not len(bw):return dict(identity_loss=1.,composition_tv=1.,empty_both=False)
    d=cdist(ax,bx);reward=(min(len(ax),len(bx))+1)*radius;cost=np.zeros((len(ax)+len(bx),len(ax)+len(bx)));cost[:len(ax),:len(bx)]=np.where(d<=radius,d-reward,reward)
    i,j=linear_sum_assignment(cost);ok=(i<len(ax))&(j<len(bx));i,j=i[ok],j[ok];ok=d[i,j]<=radius;i,j=i[ok],j[ok];correct=int((aw[i]==bw[j]).sum());negative=np.union1d(a.get('negative_ids',[]),b.get('negative_ids',[]));ac=aw[~np.isin(aw,negative)];bc=bw[~np.isin(bw,negative)];keys=np.union1d(ac,bc);ca=np.array([(ac==w).sum() for w in keys])/max(len(ac),1);cb=np.array([(bc==w).sum() for w in keys])/max(len(bc),1)
    return dict(identity_loss=1-2*correct/(len(aw)+len(bw)),composition_tv=float(np.abs(ca-cb).sum()/2) if len(ac) and len(bc) else (0. if not len(ac) and not len(bc) else 1.),composition_undefined=not len(ac) or not len(bc),empty_both=False)


def output_risk(catalogs, compatibility_result, incumbent, conventional_indices, min_calls=10, identity_threshold=.05, composition_threshold=.05):
    if compatibility_result['abstained']:return dict(abstained=True,reason=compatibility_result['reason'])
    base=catalogs[incumbent]
    if np.sum(~np.isin(base['words'],base.get('negative_ids',[])))<min_calls:return dict(abstained=True,reason='Insufficient incumbent calls')
    plausible=np.flatnonzero(compatibility_result['compatible']).tolist();rows=[dict(index=k,**catalog_distance(base,catalogs[k])) for k in plausible]
    identity=max(r['identity_loss'] for r in rows);tv=max(r['composition_tv'] for r in rows)
    conventional=[catalog_distance(base,catalogs[k]) for k in conventional_indices]
    ordinary=max([r['identity_loss'] for r in conventional],default=0.)
    all_rows=[catalog_distance(base,c) for c in catalogs]
    # If every searched candidate yields the same catalog, this still cannot
    # guarantee RNA correctness or correctness outside the transform bank.
    return dict(abstained=False,identity_risk=identity,composition_risk=tv,flagged=identity>=identity_threshold or tv>composition_threshold,compatible_candidates=len(plausible),incumbent_compatible=bool(compatibility_result['compatible'][incumbent]),ordinary_multistart_identity_disagreement=ordinary,ordinary_multistart_composition_disagreement=max([r['composition_tv'] for r in conventional],default=0.),unfiltered_bank_identity_disagreement=max(r['identity_loss'] for r in all_rows),unfiltered_bank_composition_disagreement=max(r['composition_tv'] for r in all_rows),comparisons=rows,molecular_correctness_certified=False)
