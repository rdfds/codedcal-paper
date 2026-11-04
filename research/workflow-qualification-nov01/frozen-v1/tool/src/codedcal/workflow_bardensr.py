"""Native BarDensr workflow adapter with explicit conditioned-input scope.

The optional dependency is imported only by this adapter. No policy receives
injected maps or native-reference labels. All map gauges are fixed before decode.
"""
from __future__ import annotations
import time
import hashlib
import copy
import numpy as np
from scipy.ndimage import maximum_filter
from scipy.signal import fftconvolve
from .workflow_validation import Catalog,WorkflowResult,apply_corrections,canonical,warp


def bounded_translation(fixed,moving,maximum):
    fixed=np.asarray(fixed,float);moving=np.asarray(moving,float)
    corr=fftconvolve(fixed,moving[tuple(slice(None,None,-1) for _ in moving.shape)],mode='full')
    middle=np.array(moving.shape)-1;limit=np.minimum(np.array(maximum,int),np.array(moving.shape)-1)
    slices=tuple(slice(int(c-m),int(c+m+1)) for c,m in zip(middle,limit))
    local=corr[slices];p=np.array(np.unravel_index(np.argmax(local),local.shape));result=(p-limit).astype(float)
    for axis in range(len(result)):
        if 0<p[axis]<local.shape[axis]-1:
            q=p.copy();q[axis]-=1;l=float(local[tuple(q)]);q[axis]+=2;r=float(local[tuple(q)]);c=float(local[tuple(p)])
            denominator=l-2*c+r
            if denominator<0:result[axis]+=float(np.clip(.5*(l-r)/denominator,-.5,.5))
    return result


class BardensrAdapter:
    def __init__(self,config=None):
        self.config=dict(config or {});self.decoder=self.config.get('decoder','correlation')
        if self.decoder not in ('correlation','regression'):raise ValueError('Unknown native decoder')
        import bardensr
        self.native=bardensr
        self.cache={};self.geometry_cache={}

    def context(self,case):
        z,y,x=case.images.shape[1:]
        return tuple(slice(*r) for r in case.metadata.get('context',[[0,z],[16,y-16],[16,x-16]]))

    def score_region(self,case):
        z,y,x=case.images.shape[1:]
        return tuple(tuple(r) for r in case.metadata.get('scoring_region',[[0,z],[24,y-24],[24,x-24]]))

    def objective(self,images,book):
        # Actual image/codebook residual diagnostic, not a noise likelihood.
        a=images.reshape(len(images),-1).astype(np.float32);b=book/np.sqrt(np.square(book).sum(axis=0))
        best=np.zeros(a.shape[1],np.float32)
        for j in range(0,len(b.T),128):best=np.maximum(best,(b[:,j:j+128].T@a).max(axis=0))
        return max(0.,1-float(np.square(best.astype(float)).sum())/max(float(np.square(a.astype(float)).sum()),1e-20))

    def group_fit(self,case,common):
        """Generalized retained group search, checked against the native loss.

        Two fixed initializations and the same fixed integer grid per color split.
        This is an existing comparator, not the proposed qualification method.
        """
        import itertools
        from bardensr.registration import lowrankregistration as LR
        region=tuple(slice(*r) for r in case.metadata.get('calibration_region',[[0,case.images.shape[1]],[0,case.images.shape[2]],[0,case.images.shape[3]]]))
        images=case.images[(slice(None),)+region];alpha=case.support[(slice(None),)+region]
        images=np.where(alpha>=1-1e-6,images,0).max(axis=1).astype(np.float32)
        colors=np.unique(case.colors)
        if len(colors) not in (2,4):raise ValueError('Group-search comparator qualified for two/four colors only')
        partitions=[(colors[0],)] if len(colors)==2 else [(colors[0],c) for c in colors[1:]]
        normalized=(case.codebook/np.sqrt(np.square(case.codebook).sum(axis=0))).astype(np.float32)
        denominator=max(float(np.square(images.astype(float)).sum()),1e-20);best=None;fits=[]
        ysize,xsize=images.shape[1:];inner=(ysize+40,xsize+40)
        for start_name,initial in [('zero',np.zeros_like(common[:,1:])),('round',np.rint(common[:,1:]))]:
            initial=initial.astype(int)
            if np.max(np.abs(initial))>10:continue
            shifted=warp(np.pad(images,((0,0),(30,30),(30,30)))[:,None],np.c_[np.zeros(len(images)),initial])[:,0]
            localbest=None
            for partition in partitions:
                mask=np.isin(case.colors,partition)
                plus=(normalized[mask].T@shifted[mask].reshape(mask.sum(),-1)).reshape(-1,ysize+60,xsize+60)
                minus=(normalized[~mask].T@shifted[~mask].reshape((~mask).sum(),-1)).reshape(-1,ysize+60,xsize+60)
                buffer=np.empty((case.codebook.shape[1],)+inner,np.float32)
                def evaluate(y,x):
                    np.add(plus[:,10-y:10-y+inner[0],10-x:10-x+inner[1]],minus[:,10+y:10+y+inner[0],10+x:10+x+inner[1]],out=buffer)
                    maximum=buffer.max(axis=0).astype(float);loss=1-float(np.square(maximum).sum())/denominator
                    return (loss,y*y+x*x,tuple(partition),int(y),int(x))
                found=min(evaluate(y,x) for y in range(-8,9,2) for x in range(-8,9,2))
                y,x=found[-2:];found=min([found]+[evaluate(y+dy,x+dx) for dy in [-1,0,1] for dx in [-1,0,1]])
                if localbest is None or found<localbest:localbest=found
            loss,_,partition,y,x=localbest;sign=np.where(np.isin(case.colors,partition),1,-1)
            maps=initial+sign[:,None]*np.array([y,x])
            actual=float(LR.calc_loss(images.astype(float),case.codebook,-20-maps,zero_padding=20))
            if abs(actual-loss)>=2e-6:raise ValueError('Native group-loss equivalence failed')
            fit=(loss,float(np.square(maps).sum()),start_name,maps,actual)
            fits.append(dict(start=start_name,projected_loss=loss,native_loss=actual))
            if best is None or fit[:3]<best[:3]:best=fit
        if best is None:raise ValueError('No bounded group-search start')
        return np.c_[np.zeros(len(images)),best[3]],dict(calibration_fit_source='two_start_balanced_color_group_search',
            integer_only_policy=True,group_loss=best[4],native_loss_equivalence_checked=True,group_families=fits)

    def calibrate(self,case,policy):
        f,z,y,x=case.images.shape;maximum=(min(2,z-1),8,8)
        if policy=='identity':return np.zeros((f,3)),{'calibration_fit_source':'identity_control'}
        region=tuple(slice(*r) for r in case.metadata.get('calibration_region',[[0,z],[0,y],[0,x]]))
        frames=case.images[(slice(None),)+region].astype(np.float32)
        def round_maps():
            groups=[np.flatnonzero(case.rounds==r) for r in np.unique(case.rounds)]
            merged=np.stack([frames[g].sum(axis=0) for g in groups]);ref=merged[0]
            out=np.zeros((f,3))
            for g,image in zip(groups,merged):out[g]=bounded_translation(ref,image,maximum)
            return canonical(out)
        common=round_maps()
        if policy=='round':return common,{'calibration_fit_source':'bounded_round_image_correlation'}
        if policy=='group_search':
            if not hasattr(self,'group_fit'):raise ValueError('This adapter has no qualified group-search component')
            maps,info=self.group_fit(case,common)
            return canonical(maps),info
        initial=[('zero',np.zeros_like(common)),('round',common)]
        ref=frames[0];first=np.stack([bounded_translation(ref,a,maximum) for a in frames]);initial.append(('first_readout',canonical(first)))
        candidates=[]
        for label,maps in initial:
            m=maps.copy()
            for _ in range(5):
                aligned=warp(frames,m);template=aligned.mean(axis=0)
                delta=np.stack([bounded_translation(template,a,maximum) for a in aligned]);new=canonical(m+delta)
                if np.max(np.abs(new[:,0]))>4 or np.max(np.abs(new[:,1:]))>12:break
                if np.max(np.abs(new-m))<.05:break
                m=new
            aligned=warp(frames,m);loss=float(np.square(aligned-aligned.mean(axis=0)).sum(dtype=np.float64)/max(np.square(aligned).sum(dtype=np.float64),1e-20))
            candidates.append((loss,float(np.square(m).sum()),label,m))
        image_best=min(candidates,key=lambda q:q[:3])
        if policy=='frame_template':return image_best[3],dict(calibration_fit_source='three_start_full_volume_template',selected_start=image_best[2],template_loss=image_best[0])
        if policy=='barcode_model':
            context=self.context(case);data=case.images[(slice(None),)+context].astype(float);models=[];failures=[]
            projected=bool(case.metadata.get('model_projection',False))
            if projected:data=data.max(axis=1,keepdims=True)
            for label,guess in [('zero',np.zeros_like(common)),('round',common),('template',image_best[3])]:
                try:
                    guess=guess.copy()
                    if projected:guess[:,0]=0
                    pad=10 if data.shape[1]==1 else (3,10,10)
                    fit=self.native.registration.find_translations_using_model(data,case.codebook,maximum_wiggle=pad,niter=50,initial_guess=-guess)
                    maps=canonical(-np.asarray(fit.corrections));loss=float(np.asarray(fit.losses)[-1])
                    if not np.isfinite(maps).all() or not np.isfinite(loss):raise ValueError('Nonfinite author fit')
                    models.append((loss,float(np.square(maps).sum()),label,maps))
                except Exception as e:failures.append(dict(start=label,error=type(e).__name__+': '+str(e)))
            if not models:raise RuntimeError('All native barcode-model starts failed: '+str(failures))
            best=min(models,key=lambda q:q[:3]);return best[3],dict(calibration_fit_source='actual_native_three_start_model',selected_start=best[2],model_loss=best[0],failed_starts=failures,fit_maxZ_projection=projected)
        raise ValueError('Unknown calibration policy: '+policy)

    def decode(self,case):
        context=self.context(case);a=case.images[(slice(None),)+context]
        if self.decoder=='correlation':density=self.native.spot_calling.estimate_density_singleshot(a,case.codebook,.05)
        else:density,_=self.native.spot_calling.estimate_density_iterative(a.astype(float),case.codebook,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
        density=np.asarray(density,np.float32)
        cut=np.float32(.775 if self.decoder=='correlation' else .20*float(density.max()))
        peaks=(density==maximum_filter(density,size=(3,3,3,1),mode='constant',cval=np.finfo(np.float32).min))&(density>cut)
        z,y,x,word=np.where(peaks);coords=np.c_[z,y,x]+np.array([s.start or 0 for s in context])
        region=self.score_region(case);keep=np.all((coords>=np.array([r[0] for r in region]))&(coords<np.array([r[1] for r in region])),axis=1)
        coords=coords[keep];word=word[keep];order=np.lexsort((coords[:,2],coords[:,1],coords[:,0],word))
        catalog=Catalog(coords[order],case.class_ids[word[order]])
        return catalog,dict(native_kernel='BarDensr '+self.decoder,absolute_cut=float(cut),peak_ties_preserved=True)

    def run(self,case,policy):
        start=time.perf_counter()
        h=hashlib.sha256()
        for a in [case.images,case.support,case.codebook,case.rounds,case.colors,case.class_ids,case.blank_classes]:
            a=np.ascontiguousarray(a);h.update(str((a.dtype.str,a.shape)).encode());h.update(memoryview(a).cast('B'))
        import json
        h.update(json.dumps({k:v for k,v in self.config.items() if k!='decoder'},sort_keys=True).encode())
        h.update(json.dumps(case.metadata,sort_keys=True).encode());h.update(policy.encode());key=h.hexdigest()
        if key in self.cache:
            result=copy.deepcopy(self.cache[key])
            for k in ['actual_calibration_executions','actual_decoder_executions','actual_workflow_executions']:result.diagnostics[k]=0
            for k in ['calibration_seconds','decoder_seconds']:result.diagnostics[k]=0.
            result.diagnostics.update(exact_result_cache_hit=True,wall_seconds=time.perf_counter()-start)
            return result
        fit_start=time.perf_counter()
        reused=key in self.geometry_cache
        if reused:maps,fit=copy.deepcopy(self.geometry_cache[key])
        else:
            maps,fit=self.calibrate(case,policy);self.geometry_cache[key]=copy.deepcopy((maps,fit))
        fit_seconds=0. if reused else time.perf_counter()-fit_start
        out,maps=apply_corrections(case,maps);context=self.context(case)
        coverage_region=tuple(slice(*r) for r in case.metadata.get('coverage_region',[[s.start,s.stop] for s in context]))
        coverage=float(np.all(out.support[(slice(None),)+coverage_region]>=1-1e-6,axis=0).mean())
        begin=time.perf_counter();catalog,decode_info=self.decode(out);decode_seconds=time.perf_counter()-begin
        voxel_coverage=coverage
        if 'candidate_support_coverage' in decode_info:coverage=float(decode_info['candidate_support_coverage'])
        nc=int(np.sum(~np.isin(catalog.classes,case.blank_classes)));nb=len(catalog.classes)-nc
        codes=np.unique(case.class_ids);b=len(np.intersect1d(codes,case.blank_classes));proxy=nb/max(nc,1)*(len(codes)-b)/b if b else None
        # Downsample only this diagnostic, never the actual decoder input.
        diagnostic=out.images[(slice(None),)+context][:,::max(1,out.images.shape[1]//6),::2,::2]
        result=WorkflowResult(catalog,maps,coverage,dict(**fit,**decode_info,calibration_seconds=fit_seconds,decoder_seconds=decode_seconds,
            fit_loss=self.objective(diagnostic,case.codebook),blank_proxy=proxy,nonblank_calls=nc,blank_calls=nb,
            voxel_support_coverage=voxel_coverage,support_measure='source_candidate_eligibility' if 'candidate_support_coverage' in decode_info else 'all_readout_scoring_voxels',
            actual_calibration_executions=int(not reused),actual_decoder_executions=1,actual_workflow_executions=1,
            geometry_cache_hit=reused,exact_result_cache_hit=False,wall_seconds=time.perf_counter()-start))
        self.cache[key]=copy.deepcopy(result)
        return result


def make_adapter(config=None):return BardensrAdapter(config)
