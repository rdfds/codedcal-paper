"""New local-flow and normalized object-kernel qualification; no native scoring."""
import importlib.metadata as metadata
import numpy as np
from scipy.ndimage import gaussian_filter
import adapter as A
import local_trace as L
import volume_native as V

def main():
    A.guard();codes,negative,names,words=A.book();rng=np.random.default_rng(61008)
    a=np.zeros((30,128,128),np.float32);a[14,64,64]=1
    field=np.zeros(a.shape+(3,),np.float64);field[...,0]=2;field[...,1]=-3;field[...,2]=1
    assert np.array_equal(L.warp(a,field),A.integer_warp(a,(-1,3,-2)))
    flow,diag=L.fit(a,a);assert np.count_nonzero(flow)==0 and diag['identity']
    genes=rng.choice(998,60,replace=False);centers=np.array([(y,x) for y in range(29,100,10) for x in range(29,100,10)])[:60];depth=np.array([6,14,22])[np.arange(60)%3]
    images=np.zeros((24,30,128,128),np.float32)
    for gene,(y,x),z in zip(genes,centers,depth):
        for r,c in enumerate(words[gene]):images[r*4+c,z,y,x]=1
    images=gaussian_filter(images,(0,.7,1,1));images/=images.max();support=np.ones_like(images,np.uint8)
    points,miss=V.candidates(images,support);assert len(points)==60 and miss==0
    trace=L.normalized_traces(images,points);direct=[]
    for x,y,z in points:
        row=images[:,z-1:z+2,y-2:y+3,x-2:x+3].sum(axis=(1,2,3),dtype=np.float64).reshape(6,4)
        direct.append((row/(np.sqrt((row*row).sum(axis=1,keepdims=True))+1e-6)).ravel())
    assert np.allclose(trace,np.stack(direct,axis=1),rtol=1e-6,atol=1e-7)
    results=[];lookup={(float(x),float(y),float(z)):int(j) for (y,x),z,j in zip(centers,depth,genes)}
    for decoder in A.DECODERS:
        xyz,calls,runtime=L.decode(images,support,codes,decoder);correct=sum(lookup.get(tuple(p))==int(j) for p,j in zip(xyz,calls));f1=2*correct/(60+len(calls));nb=int(negative[calls].sum());r=dict(decoder=decoder,f1=f1,coding=len(calls)-nb,negative=nb,correct=correct,runtime=runtime);results.append(r)
        A.write(L.DEST/'controls.json',dict(passed=False,synthetic=results));print('LOCAL CONTROL',decoder,r,flush=True);assert f1>=.99 and nb==0
    # A nontrivial local-flow invocation must be finite with consistent output
    # shape; source global translation controls establish integer recovery.
    fixed=gaussian_filter(a,1);moving=A.integer_warp(fixed,(0,1,1));flow,diag=L.fit(fixed,moving);assert flow.shape==a.shape+(3,) and np.isfinite(flow).all()
    before=np.square(fixed-moving).mean();after=np.square(fixed-L.warp(moving,flow)).mean();assert after<before
    A.write(L.DEST/'controls.json',dict(passed=True,synthetic=results,seed=61008,pull_warp_axes_support=True,identity_exact=True,source_trace_reference=True,known_translation_image_mse_before=float(before),known_translation_image_mse_after=float(after),local_flow_finite=True,script_sha256=A.sha(__file__),adapter_sha256=A.sha(L.__file__),development=True,independent_confirmation=False))
    paths=[A.HERE/'local_trace.py',A.HERE/'check_local_trace.py',A.HERE/'LOCAL_TRACE_PROTOCOL.md',A.HERE/'volume_native.py',L.DEST/'controls.json',A.HERE/'freeze.json',A.HERE/'preparation.json']
    base=__import__('json').loads((A.HERE/'freeze.json').read_text())['files'].copy();base.update({str(p.relative_to(A.ROOT)):A.sha(p) for p in paths})
    A.write(L.DEST/'freeze.json',dict(files=base,SimpleITK_version=metadata.version('SimpleITK'),development_pixels_already_exposed=True,independent_comparative_workflow_frozen=False,source_field='Position351'))
    print('PASSED normalized-trace and 3D local source-family controls',flush=True)

if __name__=='__main__':main()
