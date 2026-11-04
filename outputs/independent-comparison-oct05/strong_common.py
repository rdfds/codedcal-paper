"""Additional strong common-round comparator; original runner is unchanged."""
import run as study
import argparse
import time
import numpy as np

def common_iterative(observed, initial):
    frames=study.prep(observed)[:,0].reshape(7,4,200,200).sum(axis=1)
    current=initial.reshape(7,4,2).mean(axis=1)
    steps=[]
    for _ in range(10):
        template=study.warp(frames,current).mean(axis=0)
        current=np.array([study.bounded_correlation(template,frame,10) for frame in frames])
        current-=current.mean(axis=0)
        steps.append(current.tolist())
    return np.repeat(current,4,axis=0),steps

def execute(split):
    x,book,xy,words,meta=study.load(split)
    settings=study.json.loads((study.HERE/'settings-freeze.json').read_text())
    dest=study.HERE/(split+'-strong-common');dest.mkdir(exist_ok=True)
    for source in sorted((study.HERE/split).glob('*-d*-s*.json')):
        original=study.json.loads(source.read_text())
        name=original['condition'];path=dest/(name+'.json')
        if path.exists():continue
        imposed=np.array(original['imposed_yx'])
        observed=x.copy() if name.startswith('zero') else study.warp(x,imposed)
        start=time.monotonic()
        correction,steps=common_iterative(observed,np.array(original['corrections_yx']['common_round']))
        calibration_s=time.monotonic()-start+original['calibration_s']['common_round']
        start=time.monotonic();aligned=study.warp(observed,correction);warp_s=time.monotonic()-start
        result=dict(condition=name,source_sha256=meta['sha256'],
            original_condition_sha256=study.digest(source),
            amendment_sha256=study.digest(study.HERE/'STRONG_COMMON_AMENDMENT.md'),
            strong_script_sha256=study.digest(__file__),
            settings_sha256=study.digest(study.HERE/'settings-freeze.json'),
            correction_yx=correction.tolist(),steps=steps,arms={},failures={})
        for decoder in study.GRIDS:
            try:
                density,runtime=study.decode(aligned,book,decoder)
                threshold=settings['thresholds'][decoder]
                pred,ids=study.calls(density,decoder,threshold)
                key='common_iterative__'+decoder
                np.savez_compressed(dest/(name+'__'+key+'.npz'),density=density,xy=pred,words=ids)
                metrics={str(r):study.score(xy,words,pred,ids,book.shape[1],r) for r in [1.5,3,5]}
                sensitivity={}
                for multiplier in [.5,2.]:
                    pp,jj=study.calls(density,decoder,threshold*multiplier)
                    sensitivity[str(multiplier)]=study.score(xy,words,pp,jj,book.shape[1])
                result['arms'][key]=dict(metrics=metrics,threshold_sensitivity=sensitivity,
                    residual_imposed_rms=float(np.sqrt(np.mean(np.sum((correction+imposed)**2,axis=1)))),
                    runtime=dict(**runtime,calibration_s=calibration_s,warp_s=warp_s),threshold=threshold)
                print(split,name,key,'F1',round(metrics['5']['identity_f1'],4),flush=True)
            except Exception:
                result['failures'][decoder]=study.traceback.format_exc()
        study.write(path,result)
    study.write(dest/'completion.json',dict(completed_conditions=len(list(dest.glob('*-d*-s*.json')))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('split',choices=['development','validation']);a=p.parse_args()
    execute(a.split)
