"""Post-hoc development sampling repair; selection uses first-round peaks only."""
import argparse
import gc
import json
from pathlib import Path
import time
import traceback
import numpy as np
from scipy.ndimage import label
from skimage.morphology import local_maxima
import adapter as A

DEST=A.HERE/'active-plane'

def select(images,support):
    assert images.shape==(24,30,128,128) and support.shape==images.shape
    evidence=np.zeros(30,int)
    for c in range(4):
        a=images[c]
        peaks=local_maxima(a,connectivity=3,allow_borders=True)&(a>.2)
        # Count connected plateaus, not their voxel area. Lower z wins ties.
        labels,n=label(peaks,structure=np.ones((3,3,3)))
        from scipy.ndimage import center_of_mass
        for center in center_of_mass(peaks,labels,np.arange(1,n+1)):
            z,y,x=np.floor(np.asarray(center)+.5).astype(int)
            if 24<=y<104 and 24<=x<104:evidence[z]+=1
    safe=np.array([z for z in range(4,26) if np.all(support[:,z-1:z+2,24:104,24:104])])
    assert len(safe),'No supported candidate plane'
    selected=int(min(safe,key=lambda z:(-int(evidence[z]),int(z))))
    return selected,dict(selected_z=selected,reference_spot_counts_by_z=evidence.tolist(),safe_planes=safe.tolist(),selected_reference_spots=int(evidence[selected]),selection_uses_round_one_only=True,selection_uses_gene_identity_or_method_outcomes=False)

def center(images,z):
    return np.stack([A.integer_warp(frame,(14-z,0,0)) for frame in images])

def guard():
    spec=json.loads((DEST/'freeze.json').read_text())
    for rel,digest in spec['files'].items():assert A.sha(A.ROOT/rel)==digest,rel
    assert json.loads((DEST/'controls.json').read_text())['passed']
    return spec

def run(decoders=None):
    guard();codes,negative,*_=A.book();records=[];selectors=[]
    selected_decoders=A.DECODERS if decoders is None else decoders
    for k in range(5):
        path=A.CACHE/'conditioned-Position351'/f'crop-{k}.npz'
        with np.load(path) as data:images=data['images'];support=data['support']
        z,diagnostic=select(images,support);selectors.append(diagnostic|dict(crop=k))
        images=center(images,z);support=center(support,z)
        for decoder in selected_decoders:
            target=DEST/'native'/f'crop-{k}__{decoder}.json'
            if target.exists():
                result=json.loads(target.read_text());assert result['freeze_sha256']==A.sha(DEST/'freeze.json')
            else:
                try:
                    start=time.monotonic();xy,calls,runtime=A.decode(images,codes,decoder)
                    result=dict(status='completed',crop=k,decoder=decoder,selection=diagnostic,xy=xy.tolist(),words=calls.tolist(),coding=int((~negative[calls]).sum()),negative=int(negative[calls].sum()),support=float(np.all(support[:,13:16,24:104,24:104],axis=(0,1)).mean()),runtime=runtime,elapsed_s=time.monotonic()-start)
                except Exception:result=dict(status='failed',crop=k,decoder=decoder,traceback=traceback.format_exc())
                result.update(input_sha256=A.sha(path),freeze_sha256=A.sha(DEST/'freeze.json'));A.write(target,result)
            records.append(result);print('ACTIVE NATIVE',k,z,decoder,result['status'],result.get('coding'),result.get('negative'),flush=True);gc.collect()
        A.write(DEST/'selection.json',dict(field='Position351',development=True,selectors=selectors,independent_confirmation=False))
    summarize()

def summarize():
    groups={};outcomes=[]
    for decoder in A.DECODERS:
        rows=[]
        for k in range(5):
            path=DEST/'native'/f'crop-{k}__{decoder}.json'
            r=json.loads(path.read_text()) if path.exists() else dict(status='not_run',crop=k,decoder=decoder)
            rows.append(r);outcomes.append(r)
        n=sum(r.get('coding',0) for r in rows);b=sum(r.get('negative',0) for r in rows)
        complete=all(r['status']=='completed' for r in rows);proxy=b/n*998/26 if n else None;support=min(r.get('support',0) for r in rows)
        groups[decoder]=dict(completed=complete,coding=n,negative=b,scaled_proxy=proxy,support=support,qualified=complete and n>=50 and proxy<=.1 and support>=.99)
    A.write(DEST/'state.json',dict(stage='native_qualified' if all(r['qualified'] for r in groups.values()) else 'partial_or_native_gate_failure',planned=15,completed=sum(r['status']=='completed' for r in outcomes),failed=sum(r['status']=='failed' for r in outcomes),not_run=sum(r['status']=='not_run' for r in outcomes),groups=groups,development=True,adapter_qualified=all(r['qualified'] for r in groups.values()),independent_comparison_completed=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--decoder',choices=A.DECODERS);p.add_argument('--summarize',action='store_true');a=p.parse_args()
    if a.summarize:summarize()
    else:run([a.decoder] if a.decoder else None)
