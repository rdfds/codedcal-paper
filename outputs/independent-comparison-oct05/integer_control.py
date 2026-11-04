"""Post-hoc integer-only resampling sensitivity, all specified arms retained."""
import argparse
import time
import traceback
import numpy as np
import run as study
from skimage.registration import phase_cross_correlation

def imposed(kind,amplitude):
    if kind=='zero':return np.zeros((28,2),int)
    rng=np.random.default_rng(504)
    pairs=rng.integers(-amplitude,amplitude+1,size=(3,2))
    rounds=np.r_[pairs,-pairs,np.zeros((1,2),int)];rng.shuffle(rounds)
    result=np.repeat(rounds,4,axis=0)
    if kind=='color':
        c=rng.integers(-amplitude,amplitude+1,size=(2,2));colors=np.r_[c,-c];rng.shuffle(colors)
        result+=np.tile(colors,(7,1))
    assert np.array_equal(result.mean(axis=0),[0,0]);return result

def refine(images,initial):
    current=initial.copy()
    for _ in range(10):
        template=study.warp(images,current).mean(axis=0)
        current=np.array([study.bounded_correlation(template,a,10) for a in images]);current-=current.mean(axis=0)
    return current

def execute(split):
    raw,book,xy,words,source=study.load(split);x=study.prep(raw)[:,0]
    dest=study.HERE/'integer-control'/split;dest.mkdir(parents=True,exist_ok=True)
    settings=study.HERE/'settings-freeze.json';import json
    thresholds=json.loads(settings.read_text())['thresholds']
    for kind,amplitude in [('zero',0),('round',2),('round',4),('color',2),('color',4)]:
        name=f'{kind}-a{amplitude}';path=dest/(name+'.json')
        if path.exists():continue
        offsets=imposed(kind,amplitude);observed=study.warp(x,offsets)
        oracle=study.warp(observed,-offsets)
        assert np.array_equal(x[:,24:176,24:176],oracle[:,24:176,24:176])
        out=dict(condition=name,split=split,seed=504,source=source,integer_control_sha256=study.digest(__file__),
            specification_sha256=study.digest(study.HERE/'INTEGER_CONTROL.md'),original_runner_sha256=study.digest(study.HERE/'run.py'),
            offsets_yx=offsets.tolist(),added_rms=float(np.sqrt(np.mean(np.sum(offsets**2,axis=1)))),
            oracle_interior_exact=True,failures={},arms={})
        start=time.monotonic();frames=observed
        summed=frames.reshape(7,4,200,200).sum(axis=1)
        phase=np.array([phase_cross_correlation(summed[0],im,upsample_factor=50)[0] for im in summed]);phase-=phase.mean(axis=0);phase=np.repeat(phase,4,axis=0)
        common=np.repeat(refine(summed,phase[::4]),4,axis=0)
        perframe=refine(frames,phase)
        maps=dict(none=np.zeros((28,2)),common_round=phase,common_iterative=common,population_template=perframe,oracle=-offsets)
        feature=float(np.sqrt(np.mean(np.sum((perframe-common)**2,axis=1))))
        out.update(feature=feature,recommend='population_template' if feature>=2 else 'common_iterative',template_calibration_s=time.monotonic()-start)
        models={};losses={}
        for initial in ['common_round','population_template']:
            try:
                guess=np.c_[np.zeros(28),-maps[initial]].astype(float)
                fit=study.bardensr.registration.find_translations_using_model(observed[:,None].astype(float),book,maximum_wiggle=10,niter=50,initial_guess=guess)
                models[initial]=-np.asarray(fit.corrections)[:,1:];losses[initial]=np.asarray(fit.losses).tolist()
                assert np.isfinite(models[initial]).all() and np.isfinite(losses[initial]).all()
            except Exception:out['failures']['barcode_model_initial_'+initial]=traceback.format_exc()
        if models:
            selected=min(models,key=lambda k:losses[k][-1]);maps['barcode_model']=models[selected]
            out.update(model_losses=losses,model_selected_initial=selected)
        out['unrounded_corrections_yx']={k:v.tolist() for k,v in maps.items()}
        for policy in ['none','common_round','common_iterative','population_template','barcode_model','oracle']:
            if policy not in maps:continue
            correction=np.rint(maps[policy]).astype(int);aligned=study.warp(observed,correction)
            for decoder in ['correlation','regression']:
                key=policy+'__'+decoder
                try:
                    start=time.monotonic();data=aligned[:,None,16:184,16:184]
                    if decoder=='correlation':density=study.bardensr.spot_calling.estimate_density_singleshot(data,book,.01)
                    else:density,_=study.bardensr.spot_calling.estimate_density_iterative(data.astype(float),book,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
                    density=np.asarray(density,np.float32);assert np.isfinite(density).all()
                    runtime=dict(decoder_s=time.monotonic()-start)
                    pred,ids=study.calls(density,decoder,thresholds[decoder]);pred+=16
                    metrics={str(r):study.score(xy,words,pred,ids,book.shape[1],r) for r in [1.5,3,5]}
                    out['arms'][key]=dict(metrics=metrics,correction_yx=correction.tolist(),runtime=runtime,
                        residual_added_rms=float(np.sqrt(np.mean(np.sum((offsets+correction)**2,axis=1)))),
                        xy=pred.tolist(),words=ids.tolist(),density_max=float(density.max()))
                    print(split,name,key,metrics['5']['identity_f1'],flush=True)
                except Exception:out['failures'][key]=traceback.format_exc()
                study.write(path,out)
        study.write(path,out)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('split',choices=['development','validation']);a=p.parse_args();execute(a.split)
