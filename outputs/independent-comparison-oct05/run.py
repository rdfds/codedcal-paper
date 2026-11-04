"""Frozen, annotated public BARseq comparison; no molecular truth claims."""
import os
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '2')
os.environ.setdefault('TF_NUM_INTRAOP_THREADS', '4')
os.environ.setdefault('TF_NUM_INTEROP_THREADS', '1')
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback
import h5py
import numpy as np
from scipy.ndimage import shift
from scipy.signal import fftconvolve
from skimage.registration import phase_cross_correlation

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
UPSTREAM = ROOT.parent / 'bardensr-upstream.nosync'
sys.path[:0] = [str(UPSTREAM), str(ROOT / 'tool/src')]
import bardensr
from codedcal.evaluation import spatial_match

FILES = {
    'development': ('ab701a5a-2dc3-11eb-9890-0242ac110002',
                    'd78bd29c10642c7f0ee30ae2d29445e0d0b5c44a6329846fc926f0f90751fbf6'),
    'validation': ('bc_L2_MAX_M1-L2_000_002_x1500_y600_v004',
                   '520cabe0d7ac2cb7ec869219a7a477c846f69faea3d59c38090ab6048b793b72')}
GRIDS = {'correlation': [.50,.60,.65,.70,.72,.75,.80,.85,.90],
         'regression': [.03,.05,.075,.10,.15,.20,.30,.40]}

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)

def inside(xy):
    return np.all((xy >= 24) & (xy < 176), axis=1)

def score(xy, words, pred, ids, nwords, radius=5):
    keep = inside(pred)
    pred, ids = pred[keep], ids[keep]
    pairs, diag = spatial_match(xy, pred, radius)
    correct = int(np.sum(words[pairs[:,0]] == ids[pairs[:,1]]))
    total = len(xy) + len(pred)
    a = np.bincount(words, minlength=nwords)
    b = np.bincount(ids, minlength=nwords)
    constrained = 0
    for word in range(nwords):
        pp, _ = spatial_match(xy[words == word], pred[ids == word], radius)
        constrained += len(pp)
    return dict(reference_kind='author_hand_curated', n_reference=len(xy),
                calls=len(pred), correct_labels=correct, spatial_matches=len(pairs),
                identity_f1=2*correct/total if total else 0.,
                spatial_f1=2*len(pairs)/total if total else 0.,
                label_constrained_f1=2*constrained/total if total else 0.,
                wrong_label_fraction=(len(pairs)-correct)/len(pairs) if len(pairs) else None,
                composition_tv=float(np.abs(a/len(xy)-b/len(pred)).sum()/2) if len(xy) and len(pred) else None,
                matching=diag)

def prep(x):
    x = bardensr.preprocessing.minmax(x[:,None])
    x = bardensr.preprocessing.background_subtraction(x, [0,10,10])
    return bardensr.preprocessing.minmax(x)

def warp(x, offsets):
    return np.stack([shift(frame, t, order=1, mode='constant', cval=0,
                           prefilter=False) for frame,t in zip(x,offsets)])

def decode(x, book, decoder):
    start = time.monotonic()
    xx = prep(x)
    preprocessing_s = time.monotonic() - start
    start = time.monotonic()
    if decoder == 'correlation':
        density = bardensr.spot_calling.estimate_density_singleshot(xx,book,.01)
        learned = {}
    else:
        density, learned = bardensr.spot_calling.estimate_density_iterative(
            xx.astype(np.float64),book,iterations=100,psf_radius=(0,0,0),l1_penalty=0,
            estimate_codebook_gain=True)
    density = np.asarray(density, dtype=np.float32)
    if not np.isfinite(density).all():
        raise ValueError('nonfinite density')
    return density, dict(preprocessing_s=preprocessing_s,
                         decoder_s=time.monotonic()-start,
                         frame_gains=np.asarray(learned.get('frame_gains',[])).tolist())

def calls(density, decoder, threshold):
    absolute = threshold * float(density.max()) if decoder == 'regression' else threshold
    frame = bardensr.spot_calling.find_peaks(density,absolute,poolsize=(0,1,1))
    return frame[['m2','m1']].to_numpy(float), frame['j'].to_numpy(int)

def load(split):
    name, expected = FILES[split]
    path = UPSTREAM / 'bardensr/benchmarks' / (name + '.hdf5')
    assert digest(path) == expected, 'source hash mismatch'
    with h5py.File(path) as f:
        x = f['X'][:].reshape(28,200,200).astype(np.float32)
        book = f['codebook'][:].reshape(28,-1).astype(np.float32)
        status = f['rolonies/status'][:].astype('U')
        xy = np.c_[f['rolonies/m2'][:],f['rolonies/m1'][:]].astype(float)
        words = f['rolonies/j'][:].astype(int)
    mask = (status == 'good') & inside(xy)
    meta = dict(name=name, sha256=expected, shape=list(x.shape),
                n_words=book.shape[1], n_annotations=len(xy),
                statuses={s:int(np.sum(status==s)) for s in np.unique(status)},
                good_interior=int(mask.sum()),
                good_border_excluded=int(np.sum((status=='good') & ~inside(xy))))
    return x, book, xy[mask], words[mask], meta

def phase(x):
    xx = prep(x)[:,0].reshape(7,4,200,200).sum(axis=1)
    offsets = np.array([phase_cross_correlation(xx[0],frame,upsample_factor=50)[0]
                        for frame in xx])
    offsets -= offsets.mean(axis=0)
    return np.repeat(offsets,4,axis=0)

def code_model(x, book, common):
    # Upstream translations are *input sampling coordinates*. scipy shift is
    # an output displacement, so both the initialization and result are negated.
    initial = np.c_[np.zeros(28),-common].astype(np.float64)
    result = bardensr.registration.find_translations_using_model(
        prep(x).astype(np.float64),book,maximum_wiggle=10,niter=50,
        initial_guess=initial)
    return -np.asarray(result.corrections)[:,1:], np.asarray(result.losses).tolist()

def bounded_correlation(reference, moving, maximum_wiggle=10):
    """Amplitude correlation with bounded linear lags and quadratic subpixels.

    Raw phase normalization can amplify empty/noisy frequencies when readouts
    contain different source subsets. Linear correlation retains signal energy;
    the support bound is fixed before validation, not supplied by perturbation.
    """
    ref=reference.astype(np.float64)-reference.mean()
    mov=moving.astype(np.float64)-moving.mean()
    corr=fftconvolve(ref,mov[::-1,::-1],mode='full')
    center=np.array(moving.shape)-1
    lo=center-maximum_wiggle;hi=center+maximum_wiggle+1
    patch=corr[lo[0]:hi[0],lo[1]:hi[1]]
    peak=np.array(np.unravel_index(np.argmax(patch),patch.shape))+lo
    answer=(peak-center).astype(float)
    for axis in range(2):
        a=peak.copy();a[axis]-=1
        b=peak.copy();b[axis]+=1
        if abs(answer[axis])==maximum_wiggle:continue
        left,middle,right=corr[tuple(a)],corr[tuple(peak)],corr[tuple(b)]
        denominator=left-2*middle+right
        if denominator<0:
            answer[axis]+=np.clip(.5*(left-right)/denominator,-.5,.5)
    return answer

def template_alignment(x, common):
    xx = prep(x)[:,0]
    current = common.copy()
    steps = []
    for step in range(10):
        template = warp(xx,current).mean(axis=0)
        current = np.array([bounded_correlation(template,frame,maximum_wiggle=10)
                            for frame in xx])
        current -= current.mean(axis=0)
        steps.append(current.tolist())
    return current, steps

def offsets(kind, dose, seed):
    if kind == 'zero':
        return np.zeros((28,2))
    rng = np.random.default_rng(seed)
    t = rng.normal(size=(7,2)); t -= t.mean(axis=0)
    t *= dose / np.sqrt(np.mean(np.sum(t*t,axis=1)))
    t = np.repeat(t,4,axis=0)
    if kind == 'color':
        c = rng.normal(size=(4,2)); c -= c.mean(axis=0)
        c *= dose / np.sqrt(np.mean(np.sum(c*c,axis=1)))
        t += np.tile(c,(7,1))
    return t

def controls():
    rng = np.random.default_rng(410)
    xy = rng.uniform(45,120,(60,2))
    words = np.arange(60) % 7
    perm = np.roll(words,1)
    scores = {'identity':score(xy,words,xy,words,7),
              'balanced_label_corruption':score(xy,words,xy,perm,7),
              'geometry_shift':score(xy,words,xy+20,words,7),
              'inverse_geometry':score(xy,words,xy+20-20,words,7)}
    assert scores['identity']['identity_f1'] == 1
    assert scores['balanced_label_corruption']['identity_f1'] < .5
    assert scores['balanced_label_corruption']['spatial_f1'] == 1
    assert scores['balanced_label_corruption']['composition_tv'] == 0
    assert scores['geometry_shift']['spatial_f1'] < .5
    assert scores['inverse_geometry']['identity_f1'] == 1
    write(HERE/'controls.json',scores)

def qualify():
    controls()
    x, book, xy, words, meta = load('development')
    settings = dict(source=meta, protocol_sha256=digest(HERE/'PROTOCOL.md'),
                    run_sha256=digest(__file__), thresholds={}, outcomes={}, qualified=True)
    for decoder in GRIDS:
        density, runtime = decode(x,book,decoder)
        np.savez_compressed(HERE/f'development-zero-{decoder}.npz', density=density)
        curve = []
        for threshold in GRIDS[decoder]:
            pred, ids = calls(density,decoder,threshold)
            metrics = score(xy,words,pred,ids,book.shape[1])
            curve.append(dict(threshold=threshold, **metrics))
        best = max(curve,key=lambda row:(row['identity_f1'],row['threshold']))
        settings['thresholds'][decoder] = best['threshold']
        settings['outcomes'][decoder] = dict(curve=curve, selected=best, runtime=runtime)
        passed = best['identity_f1'] >= .70 and best['calls'] >= 50 and len(xy) >=50 and best['matching']['skipped_components']==0
        settings['qualified'] &= passed
        print(decoder,'qualified',passed,'threshold',best['threshold'],
              'F1',best['identity_f1'],'calls',best['calls'],flush=True)
    write(HERE/'development-qualification.json',settings)
    # This is the only setting-selection operation. Future runs use this file.
    write(HERE/'settings-freeze.json',dict(thresholds=settings['thresholds'],
        qualification_sha256=digest(HERE/'development-qualification.json'),
        protocol_sha256=digest(HERE/'PROTOCOL.md'),qualified=settings['qualified']))

def run(split):
    settings = json.loads((HERE/'settings-freeze.json').read_text())
    assert settings['qualified'], 'qualification failed; validation is blocked'
    assert settings['protocol_sha256'] == digest(HERE/'PROTOCOL.md')
    assert settings['amendment_sha256'] == digest(HERE/'PREVALIDATION_AMENDMENT.md')
    x,book,xy,words,meta = load(split)
    dest = HERE / split
    dest.mkdir(exist_ok=True)
    write(dest/'source.json',meta)
    conditions = [('zero',0,0)] + [(kind,dose,seed) for kind in ['round','color']
                                 for dose in [1,2,4] for seed in [501,502]]
    for kind,dose,seed in conditions:
        name = f'{kind}-d{dose}-s{seed}'
        result_path = dest/(name+'.json')
        if result_path.exists():
            continue
        imposed = offsets(kind,dose,seed)
        observed = x.copy() if kind == 'zero' else warp(x,imposed)
        start=time.monotonic(); common = phase(observed)
        phase_s=time.monotonic()-start
        maps={'none':np.zeros((28,2)), 'common_round':common, 'oracle':-imposed}
        policy_runtimes={'none':0.,'common_round':phase_s,'oracle':0.}
        failure={};losses={};model_fits={};model_runtime=0.
        try:
            start=time.monotonic(); maps['population_template'],template_steps=template_alignment(observed,common)
            policy_runtimes['population_template']=phase_s+time.monotonic()-start
        except Exception:
            maps.pop('population_template',None);template_steps=[]
            failure['population_template']=traceback.format_exc()
        for initial_name in ['common_round','population_template']:
            if initial_name not in maps: continue
            start=time.monotonic()
            try:
                correction, fit_losses=code_model(observed,book,maps[initial_name])
                if not np.isfinite(correction).all() or not np.isfinite(fit_losses).all():
                    raise ValueError('nonfinite model fit')
                losses[initial_name]=fit_losses
                model_fits[initial_name]=correction
            except Exception:
                failure['barcode_model_initial_'+initial_name]=traceback.format_exc()
            model_runtime += time.monotonic()-start
        selected_initial=None
        if model_fits:
            selected_initial=min(model_fits,key=lambda k:losses[k][-1])
            maps['barcode_model']=model_fits[selected_initial]
            policy_runtimes['barcode_model']=model_runtime+policy_runtimes.get('population_template',phase_s)
        result=dict(condition=name,source_sha256=meta['sha256'],
                    settings_sha256=digest(HERE/'settings-freeze.json'),
                    amendment_sha256=digest(HERE/'PREVALIDATION_AMENDMENT.md'),
                    run_sha256=digest(__file__), imposed_yx=imposed.tolist(),
                    imposed_rms=float(np.sqrt(np.mean(np.sum(imposed**2,axis=1)))),
                    corrections_yx={k:v.tolist() for k,v in maps.items()},
                    calibration_s=policy_runtimes, model_losses=losses,
                    model_selected_initial=selected_initial,template_steps=template_steps,
                    failures=failure,arms={})
        for policy, correction in maps.items():
            start=time.monotonic(); aligned=warp(observed,correction)
            warp_s=time.monotonic()-start
            for decoder in GRIDS:
                key=policy+'__'+decoder
                try:
                    density,runtime=decode(aligned,book,decoder)
                    threshold=settings['thresholds'][decoder]
                    pred,ids=calls(density,decoder,threshold)
                    np.savez_compressed(dest/(name+'__'+key+'.npz'),density=density,
                                        xy=pred,words=ids)
                    metrics={str(r):score(xy,words,pred,ids,book.shape[1],r) for r in [1.5,3,5]}
                    sensitivity={}
                    for multiplier in [.5,2.]:
                        pp,jj=calls(density,decoder,threshold*multiplier)
                        sensitivity[str(multiplier)]=score(xy,words,pp,jj,book.shape[1])
                    rms=float(np.sqrt(np.mean(np.sum((correction+imposed)**2,axis=1))))
                    result['arms'][key]=dict(metrics=metrics,threshold_sensitivity=sensitivity,
                                            residual_imposed_rms=rms,runtime=dict(
                                                **runtime,calibration_s=policy_runtimes[policy],
                                                warp_s=warp_s),threshold=threshold)
                    print(split,name,key,'F1',round(metrics['5']['identity_f1'],4),
                          'spatial',round(metrics['5']['spatial_f1'],4),flush=True)
                except Exception:
                    result['failures'][key]=traceback.format_exc()
        write(result_path,result)
    write(dest/'completion.json',dict(expected_conditions=len(conditions),
        completed_conditions=sum((dest/(f'{k}-d{d}-s{s}.json')).exists() for k,d,s in conditions),
        settings_sha256=digest(HERE/'settings-freeze.json')))

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['qualify','development','validation'])
    args=parser.parse_args()
    if args.stage == 'qualify': qualify()
    else: run(args.stage)
