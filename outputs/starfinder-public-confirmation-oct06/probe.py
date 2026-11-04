"""Exposed feasibility of loss-only selection against strong calibration alternatives."""
import gc
import json
import subprocess
import time
import traceback
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import adapter as A
import active_plane as P
import local_trace as L
import quality_catalog as Q
import global_repair as G
from bardensr.registration import lowrankregistration as LR

DEST=A.HERE/'probe'
POLICIES=['none','common_iterative','population_template','author_multistart','global_group','oracle']

def shifts(kind):
    v=np.array([4,2])
    if kind=='zero':return np.zeros((24,2),int)
    if kind=='common':return np.repeat(np.r_[np.tile(v,(3,1)),np.tile(-v,(3,1))],4,axis=0)
    return np.tile(np.r_[np.tile(v,(2,1)),np.tile(-v,(2,1))],(6,1))

def warp(volume,maps):
    return np.stack([A.integer_warp(frame,(0,int(y),int(x))) for frame,(y,x) in zip(volume,maps)])

def catalog(images,support,codes,negative,decoder,cut):
    xyz,calls,score,eligible,runtime=Q.scored(images,support,codes,decoder)
    keep=eligible&(score>=cut)
    return dict(xyz=xyz[keep].tolist(),words=calls[keep].tolist(),coding=int((~negative[calls[keep]]).sum()),negative=int(negative[calls[keep]].sum()),candidate_objects=len(calls),abstained=int((~keep).sum()),runtime=runtime)

def metrics(reference,prediction,radius):
    a=np.asarray(reference['xyz'],float).reshape(-1,3);b=np.asarray(prediction['xyz'],float).reshape(-1,3)
    aw=np.asarray(reference['words'],int);bw=np.asarray(prediction['words'],int)
    if not len(a):return dict(identity_f1=None,composition_tv=None,reference_calls=0,predicted_calls=len(b),undefined_empty_reference=True)
    if not len(b):return dict(identity_f1=0.,composition_tv=None,reference_calls=len(a),predicted_calls=0)
    dist=cdist(a,b);reward=(min(len(a),len(b))+1)*radius
    cost=np.zeros((len(a)+len(b),len(a)+len(b)));cost[:len(a),:len(b)]=np.where(dist<=radius,dist-reward,reward)
    rows,cols=linear_sum_assignment(cost);inside=(rows<len(a))&(cols<len(b));rows,cols=rows[inside],cols[inside];keep=dist[rows,cols]<=radius;rows,cols=rows[keep],cols[keep]
    correct=int((aw[rows]==bw[cols]).sum());ca=np.bincount(aw,minlength=1024)[:998];cb=np.bincount(bw,minlength=1024)[:998]
    tv=float(np.abs(ca/ca.sum()-cb/cb.sum()).sum()/2) if ca.sum() and cb.sum() else None
    return dict(identity_f1=2*correct/(len(a)+len(b)),composition_tv=tv,correct_labels=correct,reference_calls=len(a),predicted_calls=len(b),ambiguous_reference=int(((dist<=radius).sum(axis=1)>1).sum()),ambiguous_prediction=int(((dist<=radius).sum(axis=0)>1).sum()))

def fits(images,support,codes):
    plane,selected=P.select(images,support);flat=images[:,plane];alpha=support[:,plane]
    start=time.monotonic();common=G.common_fit(flat,6);frame=np.rint(A.R.refine(flat,common)).astype(int)
    maps=dict(none=np.zeros((24,2),int),common_iterative=common,population_template=frame);diag=dict(plane=selected,timing=dict(template_s=time.monotonic()-start))
    start=time.monotonic();models=[]
    for name,guess in [('zero',maps['none']),('round',common),('frame',frame)]:
        try:
            model=A.R.qualify.bardensr.registration.find_translations_using_model(flat[:,None].astype(float),codes,maximum_wiggle=10,niter=50,initial_guess=np.c_[np.zeros(24),-guess].astype(float))
            corrected=-np.asarray(model.corrections)[:,1:];loss=np.asarray(model.losses)
            assert np.isfinite(corrected).all() and np.isfinite(loss).all()
            models.append(dict(name=name,maps=np.rint(corrected).astype(int),losses=loss.tolist()))
        except Exception:models.append(dict(name=name,status='failed',traceback=traceback.format_exc()))
    good=[r for r in models if 'maps' in r]
    if good:maps['author_multistart']=min(good,key=lambda r:r['losses'][-1])['maps']
    diag['author_models']=[r|dict(maps=r['maps'].tolist()) if 'maps' in r else r for r in models];diag['timing']['author_s']=time.monotonic()-start
    start=time.monotonic()
    try:maps['global_group'],diag['global']=G.fit(flat,alpha,codes)
    except Exception:diag['global']=dict(status='failed',traceback=traceback.format_exc())
    diag['timing']['global_s']=time.monotonic()-start
    return maps,diag

def main():
    L.guard();state=json.loads((Q.DEST/'state.json').read_text());assert state['stage']=='development_thresholds_found'
    controls=json.loads((DEST/'controls.json').read_text());assert controls['passed'] and controls['probe_sha256']==A.sha(__file__)
    for name in ['probe.py','PROBE_PROTOCOL.md','quality_catalog.py']:
        path=A.HERE/name;assert subprocess.check_output(['git','show','HEAD:'+str(path.relative_to(A.ROOT))],cwd=A.ROOT)==path.read_bytes()
    cuts={d:state['groups'][d]['selected']['threshold'] for d in A.DECODERS}
    codes,negative,*_=A.book();physical=codes[:,:998];inputs=[]
    for k in [0,4]:
        path=A.CACHE/'local-trace-Position351'/f'crop-{k}.npz'
        with np.load(path) as data:inputs.append((data['images'],data['support']))
    references={}
    for decoder in A.DECODERS:references[decoder]=catalog(*inputs[1],codes,negative,decoder,cuts[decoder])
    A.write(DEST/'references.json',references)
    outcomes=[];decisions=[]
    for kind in ['zero','common','color']:
        imposed=shifts(kind);cal=warp(inputs[0][0],imposed);ca=warp(inputs[0][1],imposed);obs=warp(inputs[1][0],imposed);oa=warp(inputs[1][1],imposed)
        assert np.array_equal(warp(obs,-imposed)[:,:,16:112,16:112],inputs[1][0][:,:,16:112,16:112])
        start=time.monotonic();maps,diagnostics=fits(cal,ca,physical);maps['oracle']=-imposed
        z,validation_plane=P.select(obs,oa);losses={}
        for policy in ['common_iterative','population_template','author_multistart','global_group']:
            if policy in maps:losses[policy]=float(LR.calc_loss(obs[:,z].astype(float),physical,-20-np.array(maps[policy]),zero_padding=20))
        conventional=[p for p in losses if p!='global_group'];best=min(conventional,key=lambda p:losses[p])
        proposed='global_group' if 'global_group' in losses and losses[best]-losses['global_group']>=.005 else best
        decision=dict(kind=kind,losses=losses,best_conventional=best,selected=proposed,validation_plane=validation_plane,geometry_s=time.monotonic()-start,diagnostics=diagnostics,maps={p:np.asarray(v).tolist() for p,v in maps.items()},imposed_for_evaluation_only=imposed.tolist())
        decisions.append(decision);A.write(DEST/f'{kind}-geometry.json',decision);print('PROBE FIT',kind,losses,'SELECTED',proposed,flush=True)
        for policy in POLICIES:
            for decoder in A.DECODERS:
                try:
                    assert policy in maps,'Policy fit failed';aligned=warp(obs,maps[policy]);alpha=warp(oa,maps[policy])
                    if policy=='oracle':
                        points,_=L.V.candidates(aligned,alpha);native_points,_=L.V.candidates(*inputs[1])
                        assert np.array_equal(points,native_points)
                        assert np.array_equal(L.normalized_traces(aligned,points),L.normalized_traces(inputs[1][0],native_points))
                        pred=references[decoder]|dict(reference_reuse=True,exact_candidate_and_kernel_tensor_equality=True)
                    else:pred=catalog(aligned,alpha,codes,negative,decoder,cuts[decoder])
                    scores={str(r):metrics(references[decoder],pred,r) for r in [1.5,3,5]}
                    out=dict(status='completed',kind=kind,policy=policy,decoder=decoder,prediction=pred,metrics=scores)
                except Exception:out=dict(status='failed',kind=kind,policy=policy,decoder=decoder,traceback=traceback.format_exc())
                outcomes.append(out);A.write(DEST/f'{kind}__{policy}__{decoder}.json',out);print('PROBE',kind,policy,decoder,out['status'],out.get('metrics',{}).get('3',{}).get('identity_f1'),flush=True);gc.collect()
        A.write(DEST/'state.json',dict(stage='partial' if len(outcomes)<54 else 'development_probe_complete',planned=54,completed=sum(r['status']=='completed' for r in outcomes),failed=sum(r['status']=='failed' for r in outcomes),decision_cases=decisions,quality_cuts=cuts,exposed_development=True,independent_confirmation=False,source_3D_postperturbation_comparator_pending=True))

if __name__=='__main__':main()
