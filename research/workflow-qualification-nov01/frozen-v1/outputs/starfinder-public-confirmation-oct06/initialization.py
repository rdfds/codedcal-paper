"""Exposed 2D multistart diagnosis; no coded objective in conventional fit."""
import json,time,traceback
import numpy as np
import adapter as A
import active_plane as AP
import local_trace as L
import probe as P
import quality_catalog as Q
import recovery as R
DEST=A.HERE/'initialization'

def fit(images):
    rounds=len(images)//4
    merged=images.reshape(rounds,4,128,128).max(axis=1);ref=A.spectrum(merged[0])
    common=np.stack([A.register(ref,a) for a in merged]);common=np.repeat(np.rint(common-common.mean(axis=0)).astype(int),4,axis=0)
    ref=A.spectrum(images[0]);first=np.stack([A.register(ref,a) for a in images]);first=np.rint(first-first.mean(axis=0)).astype(int)
    starts=[('zero',np.zeros((len(images),2),int)),('round',common),('first_frame',first)];choices=[]
    for name,initial in starts:
        if np.max(np.abs(initial))>10:continue
        maps=initial.copy()
        for it in range(10):
            aligned=R.warp(images,maps);ref=A.spectrum(aligned.mean(axis=0));new=maps+np.stack([A.register(ref,a) for a in aligned]);new=np.rint(new-new.mean(axis=0)).astype(int)
            if np.max(np.abs(new))>10 or np.array_equal(new,maps):break
            maps=new
        aligned=R.warp(images,maps)[:,16:112,16:112].astype(float);loss=float(np.square(aligned-aligned.mean(axis=0)).sum()/max(np.square(aligned).sum(),1e-20))
        choices.append(dict(start=name,maps=maps,loss=loss,iterations=it+1))
    best=min(choices,key=lambda r:(r['loss'],float(np.square(r['maps']).sum()),r['start']))
    return best['maps'],dict(selected_start=best['start'],starts=[dict(start=r['start'],loss=r['loss'],iterations=r['iterations'],maps=r['maps'].tolist()) for r in choices])

def controls():
    from scipy.ndimage import gaussian_filter
    a=np.zeros((128,128),np.float32);a[30,30]=1;a[64,64]=.8;a[91,91]=.5;a=gaussian_filter(a,1)
    images=np.tile(a,(24,1,1));records=[]
    for kind in ['zero','common','color']:
        observed=R.warp(images,P.shifts(kind));maps,diag=fit(observed);assert np.array_equal(R.warp(observed,maps)[:,16:112,16:112],images[:,16:112,16:112]);records.append(dict(kind=kind,exact=True,diagnostics=diag))
    A.write(DEST/'controls.json',dict(passed=True,records=records,script_sha256=A.sha(__file__),protocol_sha256=A.sha(A.HERE/'INITIALIZATION_PROTOCOL.md')))

def main():
    L.guard();assert json.loads((DEST/'controls.json').read_text())['script_sha256']==A.sha(__file__)
    codes,negative,*_=A.book();cuts={d:r['selected']['threshold'] for d,r in json.loads((Q.DEST/'state.json').read_text())['groups'].items()};references=json.loads((P.DEST/'references.json').read_text())
    inputs=[]
    for k in [0,4]:
        with np.load(A.CACHE/f'local-trace-Position351/crop-{k}.npz') as data:inputs.append((data['images'],data['support']))
    records=[]
    for kind in ['zero','common','color']:
        imposed=P.shifts(kind);obs=P.warp(inputs[1][0],imposed);oa=P.warp(inputs[1][1],imposed)
        for fit_index,(images,support) in enumerate(inputs):
            image=P.warp(images,imposed);alpha=P.warp(support,imposed);z,plane=AP.select(image,alpha);start=time.monotonic();maps,diag=fit(image[:,z]);diag.update(geometry_s=time.monotonic()-start,plane=plane)
            aligned=P.warp(obs,maps);support=P.warp(oa,maps)
            for dec in A.DECODERS:
                prediction=P.catalog(aligned,support,codes,negative,dec,cuts[dec]);row=dict(source='STARfinder',kind=kind,fit_crop=[0,4][fit_index],decoder=dec,diagnostics=diag,prediction=prediction,metrics={str(r):P.metrics(references[dec],prediction,r) for r in [1.5,3,5]},status='completed');records.append(row);A.write(DEST/f'starfinder_{kind}_{fit_index}_{dec}.json',row);print('INIT STAR',kind,fit_index,dec,row['metrics']['3']['identity_f1'],flush=True)
    book,blank,labels=R.qualify.codebook();names=list(dict.fromkeys(labels));gene_ids=np.array([names.index(name) for name in labels])
    for k in [0,1]:
        images,support=R.inputs(k);refs={d:R.decode(images,book,d) for d in R.DECODERS}
        for kind in ['zero','common','color']:
            imposed=R.offsets(kind);obs=R.warp(images,imposed);oa=R.warp(support,imposed);start=time.monotonic();maps,diag=fit(obs);diag.update(geometry_s=time.monotonic()-start);aligned=R.warp(obs,maps);alpha=R.warp(oa,maps);exact=np.array_equal(aligned[:,16:112,16:112],images[:,16:112,16:112])
            for dec in R.DECODERS:
                if exact:xy,words,runtime=refs[dec];runtime=runtime|dict(reference_reuse=True,exact_kernel_input_equality=True)
                else:xy,words,runtime=R.decode(aligned,book,dec)
                ref,refwords,_=refs[dec];row=dict(source='BARseq3',field=k,kind=kind,decoder=dec,diagnostics=diag,exact_kernel_input=exact,quality=R.quality(words,blank,float(alpha[:,24:104,24:104].mean())),metrics={str(r):R.metrics(ref,refwords,xy,words,gene_ids,names,r) for r in [1.5,3,5]},runtime=runtime,status='completed',xy=xy.tolist(),words=words.tolist());records.append(row);A.write(DEST/f'barseq3_{k}_{kind}_{dec}.json',row);print('INIT BARSEQ',k,kind,dec,row['metrics']['3']['gene_class_identity_f1'],flush=True)
    A.write(DEST/'state.json',dict(stage='complete',planned=30,completed=len(records),development=True,independent_confirmation=False))
if __name__=='__main__':
    import sys
    if '--controls' in sys.argv:controls()
    else:main()
