"""Added conventional 3D registration controls and exposed outcomes."""
import gc,json,time,traceback
import numpy as np
import adapter as A
import local_trace as L
import quality_catalog as Q
import probe as P
DEST=A.HERE/'source-comparators'
POLICIES=['source_round_3d','frame_template_3d','frame_template_local_3d']

def move(images,maps):
    return np.stack([A.integer_warp(frame,shift) for frame,shift in zip(images,maps)])

def global_maps(images,roundwise):
    if roundwise:
        merged=images.reshape(6,4,*images.shape[1:]).max(axis=1)
        ref=A.spectrum(merged[0]);maps=np.stack([A.register(ref,a) for a in merged])
        maps=np.rint(maps-maps.mean(axis=0)).astype(int)
        return np.repeat(maps,4,axis=0)
    starts=[np.zeros((24,3),int),global_maps(images,True)]
    ref=A.spectrum(images[0]);first=np.stack([A.register(ref,a) for a in images]);starts.append(np.rint(first-first.mean(axis=0)).astype(int))
    choices=[]
    def bounded(m):return np.max(np.abs(m[:,0]))<=3 and np.max(np.abs(m[:,1:]))<=10
    for initial in starts:
        if not bounded(initial):continue
        maps=initial.copy()
        for step in range(10):
            aligned=move(images,maps);ref=A.spectrum(aligned.mean(axis=0))
            correction=np.stack([A.register(ref,a) for a in aligned]);new=maps+correction
            new=np.rint(new-new.mean(axis=0)).astype(int)
            if not bounded(new) or np.array_equal(new,maps):break
            maps=new
        aligned=move(images,maps)[:,4:26,16:112,16:112].astype(np.float64)
        loss=float(np.square(aligned-aligned.mean(axis=0)).sum()/np.maximum(np.square(aligned).sum(),1e-20))
        choices.append((loss,float(np.square(maps).sum()),maps))
    assert choices,'No finite bounded template start'
    return min(choices,key=lambda t:t[:2])[2]

def apply(images,support,policy):
    start=time.monotonic();roundwise=policy=='source_round_3d'
    maps=global_maps(images,roundwise);out=move(images,maps);alpha=move(support,maps);diagnostics=[]
    if policy!='frame_template_3d':
        fixed=out[:4].max(axis=0) if roundwise else out.mean(axis=0)
        blocks=[list(range(r*4,r*4+4)) for r in range(6)] if roundwise else [[f] for f in range(24)]
        for block in blocks:
            moving=out[block].max(axis=0) if roundwise else out[block[0]]
            flow,diag=L.fit(fixed,moving)
            for f in block:
                out[f]=L.warp(out[f],flow);alpha[f]=L.warp(alpha[f],flow)>=1-1e-6
            diagnostics.append(diag)
    return out,alpha,dict(geometry_s=time.monotonic()-start,maps_zyx=maps.tolist(),local=diagnostics,fit_on_scoring_observations=True)

def main():
    L.guard();control=json.loads((DEST/'controls.json').read_text());assert control['passed'] and control['script_sha256']==A.sha(__file__)
    codes,negative,*_=A.book();state=json.loads((Q.DEST/'state.json').read_text());cuts={d:state['groups'][d]['selected']['threshold'] for d in A.DECODERS}
    with np.load(A.CACHE/'local-trace-Position351/crop-4.npz') as data:base=data['images'];support=data['support']
    references=json.loads((P.DEST/'references.json').read_text());results=[]
    for kind in ['zero','common','color']:
        imposed=P.shifts(kind);observed=P.warp(base,imposed);alpha=P.warp(support,imposed)
        for policy in POLICIES:
            start=time.monotonic()
            try:out,oa,diagnostics=apply(observed,alpha,policy);error=None
            except Exception:error=traceback.format_exc()
            for decoder in A.DECODERS:
                try:
                    assert error is None,error
                    prediction=P.catalog(out,oa,codes,negative,decoder,cuts[decoder]);scores={str(r):P.metrics(references[decoder],prediction,r) for r in [1.5,3,5]}
                    row=dict(status='completed',kind=kind,policy=policy,decoder=decoder,prediction=prediction,metrics=scores,diagnostics=diagnostics)
                except Exception:row=dict(status='failed',kind=kind,policy=policy,decoder=decoder,traceback=traceback.format_exc())
                row.update(development=True,independent_confirmation=False,controls_sha256=A.sha(DEST/'controls.json'),script_sha256=A.sha(__file__))
                A.write(DEST/f'{kind}__{policy}__{decoder}.json',row);results.append(row)
                print('SOURCE COMPARATOR',kind,policy,decoder,row['status'],row.get('metrics',{}).get('3',{}).get('identity_f1'),flush=True);gc.collect()
            A.write(DEST/'state.json',dict(stage='complete' if len(results)==27 else 'partial',planned=27,completed=sum(r['status']=='completed' for r in results),failed=sum(r['status']=='failed' for r in results),development=True,independent_confirmation=False))
if __name__=='__main__':main()
