"""Full-depth source candidates with exact dominant-color / actual numerical kernels."""
import argparse
import gc
import json
import subprocess
import time
import traceback
import numpy as np
from scipy.ndimage import center_of_mass,label
from skimage.morphology import local_maxima
import adapter as A

DEST=A.HERE/'volume'

def candidates(images,support):
    points=[];unsupported=0
    for c in range(4):
        peaks=local_maxima(images[c],connectivity=3,allow_borders=True)&(images[c]>.2)
        labeled,n=label(peaks,structure=np.ones((3,3,3)))
        for centroid in center_of_mass(peaks,labeled,np.arange(1,n+1)):
            z,y,x=np.floor(np.asarray(centroid)+.5).astype(int)
            if not (4<=z<26 and 24<=y<104 and 24<=x<104):continue
            if not np.all(support[:,z-1:z+2,y-2:y+3,x-2:x+3]):unsupported+=1;continue
            points.append((x,y,z))
    return np.asarray(points,int).reshape(-1,3),unsupported

def traces(images,points):
    values=np.empty((24,len(points)),np.float32)
    for j,(x,y,z) in enumerate(points):values[:,j]=images[:,z-1:z+2,y-2:y+3,x-2:x+3].sum(axis=(1,2,3),dtype=np.float64)/75
    return values

def decode(images,support,codes,decoder):
    start=time.monotonic();points,unsupported=candidates(images,support);values=traces(images,points)
    runtime=dict(candidates=len(points),unsupported_candidates=unsupported,candidate_trace_s=time.monotonic()-start)
    if not len(points):return np.empty((0,3)),np.empty(0,int),runtime|dict(elapsed_s=time.monotonic()-start,abstention='no supported source candidates')
    if decoder=='source_spot':
        cube=values.T.reshape(-1,6,4); maxima=cube.max(axis=2)
        unique=((cube==maxima[:,:,None]).sum(axis=2)==1).all(axis=1);words=cube.argmax(axis=2)
        lookup={tuple(w):j for j,w in enumerate(codes.reshape(6,4,-1).argmax(axis=1).T)}
        keep=np.array([i for i,w in enumerate(words) if unique[i] and tuple(w) in lookup],int)
        calls=np.array([lookup[tuple(words[i])] for i in keep],int);xyz=points[keep]
    else:
        data=values[:,None,None,:]
        if decoder=='correlation':density=A.R.qualify.bardensr.spot_calling.estimate_density_singleshot(data,codes,.05)
        else:density,_=A.R.qualify.bardensr.spot_calling.estimate_density_iterative(data.astype(float),codes,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
        density=np.asarray(density,np.float32).reshape(len(points),-1);assert np.isfinite(density).all()
        threshold=float(.775 if decoder=='correlation' else .2*density.max())
        obj,calls=np.where(density>threshold);xyz=points[obj]
        runtime.update(absolute_cut=threshold,density_max=float(density.max()),density_sha256=__import__('hashlib').sha256(density.tobytes()).hexdigest(),multiple_calls_per_object=int(len(obj)-len(np.unique(obj))))
    order=np.lexsort((xyz[:,0],xyz[:,1],xyz[:,2],calls))
    return xyz[order].astype(float),calls[order],runtime|dict(elapsed_s=time.monotonic()-start)

def guard():
    A.guard();spec=json.loads((DEST/'freeze.json').read_text())
    for rel,digest in spec['files'].items():
        assert A.sha(A.ROOT/rel)==digest,rel
        assert subprocess.check_output(['git','show','HEAD:'+rel],cwd=A.ROOT)==(A.ROOT/rel).read_bytes(),rel
    assert json.loads((DEST/'controls.json').read_text())['passed']

def main():
    guard();codes,negative,*_=A.book();records=[]
    for k in range(5):
        path=A.CACHE/'conditioned-Position351'/f'crop-{k}.npz'
        with np.load(path) as d:images=d['images'];support=d['support']
        fraction=float(np.all(support[:,4:26,24:104,24:104],axis=(0,1)).mean())
        for decoder in A.DECODERS:
            target=DEST/'native'/f'crop-{k}__{decoder}.json'
            if target.exists():r=json.loads(target.read_text());assert r['freeze_sha256']==A.sha(DEST/'freeze.json')
            else:
                try:
                    xyz,calls,runtime=decode(images,support,codes,decoder)
                    r=dict(status='completed',crop=k,decoder=decoder,xyz=xyz.tolist(),words=calls.tolist(),coding=int((~negative[calls]).sum()),negative=int(negative[calls].sum()),support=fraction,runtime=runtime)
                except Exception:r=dict(status='failed',crop=k,decoder=decoder,traceback=traceback.format_exc())
                r.update(input_sha256=A.sha(path),freeze_sha256=A.sha(DEST/'freeze.json'));A.write(target,r)
            records.append(r);print('VOLUME',k,decoder,r['status'],r.get('coding'),r.get('negative'),flush=True);gc.collect()
    groups={}
    for decoder in A.DECODERS:
        rows=[r for r in records if r['decoder']==decoder];n=sum(r.get('coding',0) for r in rows);b=sum(r.get('negative',0) for r in rows);proxy=b/n*998/26 if n else None
        complete=all(r['status']=='completed' for r in rows);s=min(r.get('support',0) for r in rows)
        groups[decoder]=dict(completed=complete,coding=n,negative=b,scaled_proxy=proxy,support=s,qualified=complete and n>=50 and proxy<=.1 and s>=.99,runtime_s=sum(r.get('runtime',{}).get('elapsed_s',0) for r in rows))
    A.write(DEST/'state.json',dict(stage='native_qualified' if all(r['qualified'] for r in groups.values()) else 'native_gate_failure',planned=15,completed=sum(r['status']=='completed' for r in records),failed=sum(r['status']=='failed' for r in records),groups=groups,adapter_qualified=all(r['qualified'] for r in groups.values()),development=True,independent_comparison_completed=False))

if __name__=='__main__':main()
