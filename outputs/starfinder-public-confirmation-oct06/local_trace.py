"""Source-informed 3D local / normalized-trace development adapter."""
import gc
import importlib.metadata as metadata
import json
import subprocess
import time
import traceback
import numpy as np
import SimpleITK as sitk
import adapter as A
import volume_native as V

sitk.ProcessObject.SetGlobalDefaultNumberOfThreads(2)
DEST=A.HERE/'local-trace'

def warp(a,field):
    image=sitk.GetImageFromArray(a.astype(np.float32));flow=sitk.GetImageFromArray(field.astype(np.float64),isVector=True)
    transform=sitk.DisplacementFieldTransform(flow)
    return sitk.GetArrayFromImage(sitk.Resample(image,image,transform,sitk.sitkLinear,0.,image.GetPixelID()))

def pyramid(image,factor):
    size=np.array(image.GetSize());newsize=np.maximum(4,np.ceil(size/factor).astype(int))
    spacing=(size-1)/(newsize-1)*np.array(image.GetSpacing())
    smooth=sitk.DiscreteGaussian(image,(factor/2)**2) if factor>1 else image
    return sitk.Resample(smooth,newsize.tolist(),sitk.Transform(3,sitk.sitkIdentity),sitk.sitkLinear,image.GetOrigin(),spacing.tolist(),image.GetDirection(),0.,image.GetPixelID())

def fit(fixed,moving):
    if np.array_equal(fixed,moving):return np.zeros(fixed.shape+(3,),np.float64),dict(identity=True)
    f=sitk.GetImageFromArray(fixed.astype(np.float32));m=sitk.GetImageFromArray(moving.astype(np.float32));flow=None;log=[]
    for factor in [8,4,2,1]:
        fi=pyramid(f,factor);mi=pyramid(m,factor)
        if flow is None:
            flow=sitk.Image(fi.GetSize(),sitk.sitkVectorFloat64,3);flow.CopyInformation(fi)
        else:flow=sitk.Resample(flow,fi,sitk.Transform(3,sitk.sitkIdentity),sitk.sitkLinear,0.,sitk.sitkVectorFloat64)
        solver=sitk.DemonsRegistrationFilter();solver.SetNumberOfIterations(10);solver.SetStandardDeviations(1.)
        flow=solver.Execute(fi,mi,flow);log.append(dict(factor=factor,iterations=solver.GetElapsedIterations(),metric=solver.GetMetric(),rms_change=solver.GetRMSChange()))
    flow=sitk.Resample(flow,f,sitk.Transform(3,sitk.sitkIdentity),sitk.sitkLinear,0.,sitk.sitkVectorFloat64)
    out=sitk.GetArrayFromImage(flow);assert np.isfinite(out).all()
    return out,dict(levels=log,maximum_displacement=float(np.linalg.norm(out,axis=-1).max()),field_sha256=__import__('hashlib').sha256(out.tobytes()).hexdigest())

def normalized_traces(images,points):
    raw=V.traces(images,points)*75
    cube=raw.T.reshape(-1,6,4).astype(np.float64)
    cube/=np.sqrt((cube*cube).sum(axis=2,keepdims=True))+1e-6
    return cube.reshape(-1,24).T.astype(np.float32)

def decode(images,support,codes,decoder):
    start=time.monotonic();points,unsupported=V.candidates(images,support);values=normalized_traces(images,points)
    runtime=dict(candidates=len(points),unsupported_candidates=unsupported,candidate_trace_s=time.monotonic()-start)
    if not len(points):return np.empty((0,3)),np.empty(0,int),runtime|dict(elapsed_s=time.monotonic()-start,abstention='no supported candidates')
    if decoder=='source_spot':
        cube=values.T.reshape(-1,6,4);maxima=cube.max(axis=2);unique=((cube==maxima[:,:,None]).sum(axis=2)==1).all(axis=1)
        lookup={tuple(w):j for j,w in enumerate(codes.reshape(6,4,-1).argmax(axis=1).T)}
        words=cube.argmax(axis=2);keep=np.array([i for i,w in enumerate(words) if unique[i] and tuple(w) in lookup],int)
        calls=np.array([lookup[tuple(words[i])] for i in keep],int);xyz=points[keep]
    else:
        data=values[:,None,None,:]
        if decoder=='correlation':density=A.R.qualify.bardensr.spot_calling.estimate_density_singleshot(data,codes,.05)
        else:density,_=A.R.qualify.bardensr.spot_calling.estimate_density_iterative(data.astype(float),codes,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
        density=np.asarray(density,np.float32).reshape(len(points),-1);assert np.isfinite(density).all()
        cut=float(.775 if decoder=='correlation' else .2*density.max());obj,calls=np.where(density>cut);xyz=points[obj]
        runtime.update(absolute_cut=cut,density_max=float(density.max()),multiple_calls_per_object=int(len(obj)-len(np.unique(obj))),density_sha256=__import__('hashlib').sha256(density.tobytes()).hexdigest())
    order=np.lexsort((xyz[:,0],xyz[:,1],xyz[:,2],calls))
    return xyz[order].astype(float),calls[order],runtime|dict(elapsed_s=time.monotonic()-start)

def guard():
    A.guard();spec=json.loads((DEST/'freeze.json').read_text())
    assert metadata.version('SimpleITK')==spec['SimpleITK_version']
    for rel,digest in spec['files'].items():
        assert A.sha(A.ROOT/rel)==digest,rel
        assert subprocess.check_output(['git','show','HEAD:'+rel],cwd=A.ROOT)==(A.ROOT/rel).read_bytes(),rel
    assert json.loads((DEST/'controls.json').read_text())['passed']

def condition(k):
    directory=A.CACHE/'local-trace-Position351';directory.mkdir(exist_ok=True)
    target=directory/f'crop-{k}.npz';marker=DEST/f'crop-{k}-local.json'
    if target.exists():
        record=json.loads(marker.read_text());assert record['freeze_sha256']==A.sha(DEST/'freeze.json') and record['output_sha256']==A.sha(target);return target
    path=A.CACHE/'conditioned-Position351'/f'crop-{k}.npz'
    with np.load(path) as data:images=data['images'];support=data['support']
    fixed=images[:4].max(axis=0);out=[];alpha=[];logs=[];start=time.monotonic()
    for r in range(6):
        field,diagnostic=fit(fixed,images[r*4:r*4+4].max(axis=0))
        for c in range(4):out.append(warp(images[r*4+c],field));alpha.append(warp(support[r*4+c],field)>=1-1e-6)
        logs.append(diagnostic);print('LOCAL PREP',k,r,diagnostic.get('maximum_displacement'),flush=True)
    np.savez_compressed(target,images=np.stack(out),support=np.stack(alpha))
    A.write(marker,dict(crop=k,rounds=logs,conditioning_s=time.monotonic()-start,input_sha256=A.sha(path),output_sha256=A.sha(target),freeze_sha256=A.sha(DEST/'freeze.json')))
    return target

def main():
    guard();codes,negative,*_=A.book();records=[]
    for k in range(5):
        path=condition(k)
        with np.load(path) as data:images=data['images'];support=data['support']
        points,rejected=V.candidates(images,support);fraction=len(points)/(len(points)+rejected) if len(points)+rejected else 1.
        for decoder in A.DECODERS:
            target=DEST/'native'/f'crop-{k}__{decoder}.json'
            if target.exists():r=json.loads(target.read_text());assert r['freeze_sha256']==A.sha(DEST/'freeze.json')
            else:
                try:
                    xyz,calls,runtime=decode(images,support,codes,decoder)
                    r=dict(status='completed',crop=k,decoder=decoder,xyz=xyz.tolist(),words=calls.tolist(),coding=int((~negative[calls]).sum()),negative=int(negative[calls].sum()),support=fraction,runtime=runtime)
                except Exception:r=dict(status='failed',crop=k,decoder=decoder,traceback=traceback.format_exc())
                r.update(input_sha256=A.sha(path),freeze_sha256=A.sha(DEST/'freeze.json'));A.write(target,r)
            records.append(r);print('LOCAL NATIVE',k,decoder,r['status'],r.get('coding'),r.get('negative'),r.get('support'),flush=True);gc.collect()
    groups={}
    for decoder in A.DECODERS:
        rows=[r for r in records if r['decoder']==decoder];n=sum(r.get('coding',0) for r in rows);b=sum(r.get('negative',0) for r in rows);proxy=b/n*998/26 if n else None;s=min(r.get('support',0) for r in rows);complete=all(r['status']=='completed' for r in rows)
        groups[decoder]=dict(completed=complete,coding=n,negative=b,scaled_proxy=proxy,support=s,qualified=complete and n>=50 and proxy<=.1 and s>=.99,runtime_s=sum(r.get('runtime',{}).get('elapsed_s',0) for r in rows))
    A.write(DEST/'state.json',dict(stage='native_qualified' if all(r['qualified'] for r in groups.values()) else 'native_gate_failure',planned=15,completed=sum(r['status']=='completed' for r in records),failed=sum(r['status']=='failed' for r in records),groups=groups,adapter_qualified=all(r['qualified'] for r in groups.values()),development=True,independent_comparison_completed=False))

if __name__=='__main__':main()
