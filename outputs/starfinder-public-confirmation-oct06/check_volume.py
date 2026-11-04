"""New full-depth candidate/trace controls with fresh actual kernel executions."""
import json
import numpy as np
from scipy.ndimage import gaussian_filter
import adapter as A
import volume_native as V

def main():
    A.guard();codes,negative,names,words=A.book();rng=np.random.default_rng(61007)
    genes=rng.choice(998,60,replace=False);centers=np.array([(y,x) for y in range(29,100,10) for x in range(29,100,10)])[:60]
    images=np.zeros((24,30,128,128),np.float32)
    depth=np.array([6,14,22])[np.arange(60)%3]
    for gene,(y,x),z in zip(genes,centers,depth):
        for r,c in enumerate(words[gene]):images[r*4+c,z,y,x]=1
    images=gaussian_filter(images,(0,.7,1,1));images/=images.max();support=np.ones_like(images,np.uint8)
    points,miss=V.candidates(images,support);assert len(points)==60 and miss==0
    trace=V.traces(images,points)
    direct=np.stack([images[:,z-1:z+2,y-2:y+3,x-2:x+3].mean(axis=(1,2,3),dtype=np.float64) for x,y,z in points],axis=1).astype(np.float32)
    assert np.array_equal(trace,direct)
    outcomes=[];lookup={(float(x),float(y),float(z)):int(j) for (y,x),z,j in zip(centers,depth,genes)}
    for decoder in A.DECODERS:
        xyz,calls,runtime=V.decode(images,support,codes,decoder);correct=sum(lookup.get(tuple(p))==int(j) for p,j in zip(xyz,calls));f1=2*correct/(60+len(calls));nb=int(negative[calls].sum())
        out=dict(decoder=decoder,f1=f1,coding=len(calls)-nb,negative=nb,correct=correct,runtime=runtime);outcomes.append(out)
        A.write(V.DEST/'controls.json',dict(passed=False,synthetic=outcomes,seed=61007,adapter_sha256=A.sha(V.__file__)))
        print('VOLUME CONTROL',decoder,out,flush=True);assert f1>=.99 and nb==0
    # Source selection keeps all depths and ignores later-round spot detection.
    altered=images.copy();altered[4:]=0
    assert np.array_equal(V.candidates(altered,support)[0],points)
    _,empty,_=V.decode(np.zeros_like(images),support,codes,'regression');assert len(empty)==0
    result=dict(passed=True,seed=61007,synthetic=outcomes,reference_candidates=60,all_three_depths_retained=True,trace_mean_independent_reference=True,later_round_candidate_invariance=True,empty_abstention=True,script_sha256=A.sha(__file__),adapter_sha256=A.sha(V.__file__),development=True,independent_confirmation=False)
    A.write(V.DEST/'controls.json',result)
    paths=[A.HERE/'volume_native.py',A.HERE/'check_volume.py',A.HERE/'VOLUME_PROTOCOL.md',V.DEST/'controls.json',A.HERE/'freeze.json',A.HERE/'preparation.json']
    files=json.loads((A.HERE/'freeze.json').read_text())['files'].copy();files.update({str(p.relative_to(A.ROOT)):A.sha(p) for p in paths})
    A.write(V.DEST/'freeze.json',dict(files=files,source_field='Position351',development_pixels_already_exposed=True,independent_comparative_workflow_frozen=False,full_depth_source_defined_candidates=True))
    print('PASSED full-depth object-trace controls; no native outcome or independent claim',flush=True)

if __name__=='__main__':main()
