"""Supported single-object catalog with control-learned ambiguity thresholds."""
import json
import time
import numpy as np
from scipy.stats import chi2
import adapter as A
import local_trace as L
import volume_native as V

DEST=A.HERE/'quality-catalog'

def scored(images,support,codes,decoder):
    start=time.monotonic();points,rejected=V.candidates(images,support);values=L.normalized_traces(images,points)
    if not len(points):return points,np.empty(0,int),np.empty(0),np.empty(0,bool),dict(candidates=0,unsupported_candidates=rejected,elapsed_s=time.monotonic()-start)
    if decoder=='source_spot':
        cube=values.T.reshape(-1,6,4);sorted_colors=np.sort(cube,axis=2)
        score=(sorted_colors[:,:,-1]-sorted_colors[:,:,-2]).min(axis=1)
        lookup={tuple(w):j for j,w in enumerate(codes.reshape(6,4,-1).argmax(axis=1).T)}
        word=cube.argmax(axis=2);calls=np.array([lookup.get(tuple(w),-1) for w in word],int)
        eligible=(calls>=0)&(score>0)
    else:
        data=values[:,None,None,:]
        if decoder=='correlation':density=A.R.qualify.bardensr.spot_calling.estimate_density_singleshot(data,codes,.05)
        else:density,_=A.R.qualify.bardensr.spot_calling.estimate_density_iterative(data.astype(float),codes,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
        density=np.asarray(density,np.float32).reshape(len(points),-1);assert np.isfinite(density).all()
        ordered=np.sort(density,axis=1);first=ordered[:,-1];second=ordered[:,-2];calls=density.argmax(axis=1)
        score=(first-second) if decoder=='correlation' else (first-second)/np.maximum(density.sum(axis=1),1e-20)
        eligible=(score>0)&(first>(.775 if decoder=='correlation' else .2*density.max()))
    return points,calls,score,eligible,dict(candidates=len(points),unsupported_candidates=rejected,duplicate_candidate_coordinates=int(len(points)-len(np.unique(points,axis=0))),elapsed_s=time.monotonic()-start)

def curve(score,calls,eligible,negative):
    rows=[]
    for threshold in np.unique(np.r_[0,score[eligible]]):
        take=eligible&(score>=threshold);labels=calls[take];n=int((~negative[labels]).sum());b=int(negative[labels].sum());proxy=b/n*998/26 if n else None
        upper=float(.5*chi2.ppf(.95,2*(b+1))*998/26/n) if n else None
        rows.append(dict(threshold=float(threshold),coding=n,negative=b,scaled_proxy=proxy,negative_count_proxy_upper95=upper,qualified=n>=50 and proxy<=.1))
    passing=[r for r in rows if r['qualified']]
    selected=min(passing,key=lambda r:(-r['coding'],r['threshold'])) if passing else None
    return rows,selected

def main():
    L.guard();assert json.loads((DEST/'controls.json').read_text())['passed']
    codes,negative,*_=A.book();records=[]
    for k in range(5):
        source=A.CACHE/'local-trace-Position351'/f'crop-{k}.npz'
        with np.load(source) as d:images=d['images'];support=d['support']
        for decoder in A.DECODERS:
            target=DEST/'native'/f'crop-{k}__{decoder}.json'
            points,calls,score,eligible,runtime=scored(images,support,codes,decoder)
            r=dict(crop=k,decoder=decoder,xyz=points.tolist(),words=calls.tolist(),ambiguity_scores=score.tolist(),eligible=eligible.tolist(),runtime=runtime,input_sha256=A.sha(source),script_sha256=A.sha(__file__))
            A.write(target,r);records.append(r);print('QUALITY',k,decoder,int(eligible.sum()),flush=True)
    groups={}
    for decoder in A.DECODERS:
        rows=[r for r in records if r['decoder']==decoder]
        allscores=np.concatenate([r['ambiguity_scores'] for r in rows]);allcalls=np.concatenate([r['words'] for r in rows]).astype(int);eligible=np.concatenate([r['eligible'] for r in rows]).astype(bool)
        full,selected=curve(allscores,allcalls,eligible,negative)
        groups[decoder]=dict(selected=selected,threshold_curve=full,total_candidates=len(allcalls),eligible_before_score=int(eligible.sum()),abstained_after_score=int(len(allcalls)-selected['coding']-selected['negative']) if selected else len(allcalls))
    A.write(DEST/'state.json',dict(stage='development_thresholds_found' if all(r['selected'] is not None for r in groups.values()) else 'development_threshold_learning_failed',groups=groups,conditional_supported_object_target=True,original_native_interface_still_failed=True,development=True,independent_confirmation=False,independent_comparative_workflow_frozen=False))

if __name__=='__main__':main()
