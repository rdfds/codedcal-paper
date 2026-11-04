"""Fixed matched-covariance recovery; algorithm-output consistency, not RNA truth."""
import os
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL','2')
os.environ.setdefault('TF_NUM_INTRAOP_THREADS','2')
os.environ.setdefault('TF_NUM_INTEROP_THREADS','1')
import argparse
import hashlib
import importlib.util
from pathlib import Path
import sys
import time
import traceback
import numpy as np
from scipy.ndimage import shift,maximum_filter
from scipy.signal import fftconvolve
from skimage.registration import phase_cross_correlation
from condition import HERE,ROOT,CACHE,sha,read,write
import qualify
sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal.evaluation import spatial_match

POLICIES=['none','common_iterative','population_template','barcode_model','oracle']
DECODERS=['correlation','regression']

def offsets(kind):
    if kind=='zero':return np.zeros((40,2),int)
    rng=np.random.default_rng(505)
    v=np.array([4,2])
    if kind=='common':
        r=np.r_[np.tile(v,(5,1)),np.tile(-v,(5,1))];rng.shuffle(r);out=np.repeat(r,4,axis=0)
    else:
        colors=np.array([v,v,-v,-v]);rng.shuffle(colors);out=np.tile(colors,(10,1))
    assert np.array_equal(out.mean(axis=0),[0,0])
    assert np.array_equal(out.T@out/40,[[16,8],[8,4]])
    return out

def warp(a,maps):return np.stack([shift(im,t,order=1,mode='constant',cval=0,prefilter=False) for im,t in zip(a,maps)])

def bounded(reference,moving):
    ref=reference.astype(float)-reference.mean();mov=moving.astype(float)-moving.mean()
    corr=fftconvolve(ref,mov[::-1,::-1],mode='full');center=np.array(mov.shape)-1
    lo=center-10;hi=center+11;window=corr[lo[0]:hi[0],lo[1]:hi[1]]
    peak=np.array(np.unravel_index(np.argmax(window),window.shape))+lo;answer=(peak-center).astype(float)
    for axis in [0,1]:
        if abs(answer[axis])==10:continue
        left=peak.copy();right=peak.copy();left[axis]-=1;right[axis]+=1
        a,b,c=corr[tuple(left)],corr[tuple(peak)],corr[tuple(right)];denom=a-2*b+c
        if denom<0:answer[axis]+=np.clip(.5*(a-c)/denom,-.5,.5)
    return answer

def refine(images,initial):
    current=initial.copy()
    for _ in range(10):
        template=warp(images,current).mean(axis=0)
        current=np.array([bounded(template,im) for im in images]);current-=current.mean(axis=0)
    return current

def inputs(index):
    arrays=[];supports=[];dest=HERE/f'field-{index}'
    for r in range(10):
        p=CACHE/f'native-repair-source-{index}-round-{r}.npz';marker=read(dest/f'round-{r}-demons.json')
        assert sha(p)==marker['output_sha256'] and marker['condition_sha256']==sha(HERE/'condition.py')
        with np.load(p) as data:arrays.extend(data['images']);supports.extend(data['support'])
    return np.stack(arrays),np.stack(supports)

def geometry(observed,book):
    start=time.monotonic();summed=observed.reshape(10,4,128,128).sum(axis=1)
    initial=np.array([phase_cross_correlation(summed[0],im,upsample_factor=50)[0] for im in summed]);initial-=initial.mean(axis=0)
    shared=time.monotonic()-start
    start=time.monotonic();common=np.repeat(refine(summed,initial),4,axis=0);common_s=time.monotonic()-start
    start=time.monotonic();frame=refine(observed,np.repeat(initial,4,axis=0));frame_s=time.monotonic()-start
    maps=dict(none=np.zeros((40,2)),common_iterative=common,population_template=frame)
    losses={};models={};failures={};start=time.monotonic()
    for name,guess in [('common_round',np.repeat(initial,4,axis=0)),('population_template',frame)]:
        try:
            result=qualify.bardensr.registration.find_translations_using_model(observed[:,None].astype(float),book,
                maximum_wiggle=10,niter=50,initial_guess=np.c_[np.zeros(40),-guess].astype(float))
            models[name]=-np.asarray(result.corrections)[:,1:];losses[name]=np.asarray(result.losses).tolist()
            assert np.isfinite(models[name]).all() and np.isfinite(losses[name]).all()
        except Exception:failures['model_initial_'+name]=traceback.format_exc()
    selected=None
    if models:
        selected=min(models,key=lambda k:losses[k][-1]);maps['barcode_model']=models[selected]
    return maps,dict(shared_phase_s=shared,common_s=common_s,frame_s=frame_s,model_s=time.monotonic()-start),dict(losses=losses,selected_initial=selected,failures=failures)

def decode(a,book,decoder):
    start=time.monotonic();data=a[:,None,16:112,16:112]
    if decoder=='correlation':density=qualify.bardensr.spot_calling.estimate_density_singleshot(data,book,.05)
    else:density,_=qualify.bardensr.spot_calling.estimate_density_iterative(data.astype(float),book,iterations=100,psf_radius=(0,0,0),l1_penalty=0,estimate_codebook_gain=True)
    density=np.asarray(density,np.float32);assert np.isfinite(density).all();kernel_s=time.monotonic()-start;start=time.monotonic()
    threshold=np.float32(.775 if decoder=='correlation' else .20*float(density.max()))
    maxima=maximum_filter(density,size=(1,3,3,1),mode='constant',cval=np.finfo(np.float32).min)
    z,y,x,j=np.where((density==maxima)&(maxima>threshold));x+=16;y+=16
    inside=(x>=24)&(x<104)&(y>=24)&(y<104);x,y,j=x[inside],y[inside],j[inside]
    order=np.lexsort((x,y,j));xy=np.c_[x[order],y[order]].astype(float);j=j[order]
    return xy,j,dict(kernel_s=kernel_s,peaks_s=time.monotonic()-start,density_max=float(density.max()),absolute_cut=float(threshold),density_sha256=hashlib.sha256(density.tobytes()).hexdigest())

def quality(words,blank,support):
    nonblank=int((~blank[words]).sum());nblank=int(blank[words].sum());proxy=nblank/max(nonblank,1)*int((~blank).sum())/int(blank.sum())
    return dict(nonblank=nonblank,blanks=nblank,proxy=proxy,support=support,qualified=nonblank>=50 and proxy<=.1 and support>=.99)

def metrics(ref,refwords,pred,words,gene_ids,names,radius):
    pairs,diag=spatial_match(ref,pred,radius);refgene=gene_ids[refwords];predgene=gene_ids[words]
    correct=int((refgene[pairs[:,0]]==predgene[pairs[:,1]]).sum());seqcorrect=int((refwords[pairs[:,0]]==words[pairs[:,1]]).sum())
    n=len(ref)+len(pred);count_a=np.bincount(refgene,minlength=len(names));count_b=np.bincount(predgene,minlength=len(names))
    nonblank=np.array(['no_gene' not in name.lower() for name in names]);a=count_a[nonblank];b=count_b[nonblank]
    constrained=0
    for j in np.union1d(refgene,predgene):
        pair,_=spatial_match(ref[refgene==j],pred[predgene==j],radius);constrained+=len(pair)
    return dict(gene_class_identity_f1=2*correct/n if n else 0.,sequence_identity_f1=2*seqcorrect/n if n else 0.,
        spatial_f1=2*len(pairs)/n if n else 0.,label_constrained_f1=2*constrained/n if n else 0.,
        composition_tv=float(np.abs(a/a.sum()-b/b.sum()).sum()/2) if a.sum() and b.sum() else None,
        reference_calls=len(ref),predicted_calls=len(pred),correct_gene_classes=correct,matching=diag)

def execute(index):
    assert all(read(HERE/f'field-{i}/qualification.json')['qualified'] for i in [0,1])
    assert read(HERE.parent/'independent-comparison-oct05/integer-control/summary.json')['arms']['completed']==120
    # Public specification and implementation must be committed before outcomes.
    import subprocess
    for name in ['RECOVERY_PROTOCOL.md','recovery.py']:
        assert subprocess.check_output(['git','show','HEAD:outputs/barseq3-native-repair-oct05/'+name],cwd=ROOT)==(HERE/name).read_bytes()
    book,blank,labels=qualify.codebook();names=list(dict.fromkeys(labels));gene_ids=np.array([names.index(name) for name in labels])
    raw,rawsupport=inputs(index);dest=HERE/'recovery'/f'field-{index}';dest.mkdir(parents=True,exist_ok=True)
    references={}
    for decoder in DECODERS:
        p=dest/f'reference-{decoder}.json'
        if p.exists():ref=read(p);assert ref['recovery_sha256']==sha(__file__)
        else:
            xy,words,runtime=decode(raw,book,decoder);support=float(np.all(rawsupport[:,24:104,24:104]>=1-1e-6,axis=0).mean())
            ref=dict(xy=xy.tolist(),words=words.tolist(),quality=quality(words,blank,support),runtime=runtime,recovery_sha256=sha(__file__),protocol_sha256=sha(HERE/'RECOVERY_PROTOCOL.md'))
            write(p,ref)
        references[decoder]=ref
    if not all(r['quality']['qualified'] for r in references.values()):
        write(dest/'completion.json',dict(stage='reference_gate_failed',references={d:r['quality'] for d,r in references.items()},planned_arms=30,completed_arms=0));return
    for kind in ['zero','common','color']:
        imposed=offsets(kind);observed=warp(raw,imposed);supportin=warp(rawsupport,imposed)
        assert np.array_equal(warp(observed,-imposed)[:,16:112,16:112],raw[:,16:112,16:112])
        gp=dest/(kind+'-geometry.json')
        if gp.exists():geo=read(gp);assert geo['recovery_sha256']==sha(__file__);maps={p:np.array(v) for p,v in geo['maps'].items()}
        else:
            maps,runtime,model=geometry(observed,book);maps['oracle']=-imposed
            feature=float(np.sqrt(np.mean(np.sum((maps['population_template']-maps['common_iterative'])**2,axis=1))))
            geo=dict(kind=kind,maps={p:m.tolist() for p,m in maps.items()},imposed=imposed.tolist(),covariance=(imposed.T@imposed/40).tolist(),
                feature=feature,recommend='population_template' if feature>=2 else 'common_iterative',runtime=runtime,model=model,recovery_sha256=sha(__file__))
            write(gp,geo)
        for policy in POLICIES:
            for decoder in DECODERS:
                path=dest/f'{kind}__{policy}__{decoder}.json'
                if path.exists():continue
                try:
                    if policy not in maps:raise ValueError('Both model initializations failed')
                    correction=np.rint(maps[policy]).astype(int);start=time.monotonic();aligned=warp(observed,correction);supportout=warp(supportin,correction)
                    support=float(np.all(supportout[:,24:104,24:104]>=1-1e-6,axis=0).mean());warp_s=time.monotonic()-start
                    if support<.99:raise ValueError(f'Geometric support {support} below fixed .99')
                    ref=references[decoder]
                    if policy=='oracle' or (kind=='zero' and policy=='none'):
                        assert np.array_equal(aligned[:,16:112,16:112],raw[:,16:112,16:112])
                        xy=np.array(ref['xy'],float).reshape(-1,2);words=np.array(ref['words'],int)
                        runtime=ref['runtime']|dict(exact_reference_kernel_reused=True,actual_execution_kernel_s=0.)
                    else:
                        xy,words,runtime=decode(aligned,book,decoder);runtime['exact_reference_kernel_reused']=False
                    if policy=='oracle':assert xy.tolist()==ref['xy'] and words.tolist()==ref['words'],'Oracle changes native calls'
                    scores={str(r):metrics(np.array(ref['xy']).reshape(-1,2),np.array(ref['words'],int),xy,words,gene_ids,names,r) for r in [1.5,3,5]}
                    out=dict(status='completed',kind=kind,policy=policy,decoder=decoder,source_index=index,xy=xy.tolist(),words=words.tolist(),metrics=scores,
                        quality=quality(words,blank,support),runtime=runtime|dict(warp_s=warp_s),applied_integer_maps=correction.tolist(),
                        residual_added_rms=float(np.sqrt(np.mean(np.sum((imposed+correction)**2,axis=1)))),geometry_sha256=sha(gp),reference_sha256=sha(dest/f'reference-{decoder}.json'),recovery_sha256=sha(__file__))
                except Exception:out=dict(status='failed',kind=kind,policy=policy,decoder=decoder,source_index=index,traceback=traceback.format_exc(),recovery_sha256=sha(__file__))
                write(path,out);print(index,kind,policy,decoder,out['status'],out.get('metrics',{}).get('5',{}).get('gene_class_identity_f1'),flush=True)
    arms=[read(p) for p in dest.glob('*__*.json')]
    assert len(arms)==30
    write(dest/'completion.json',dict(stage='complete_with_retained_failures',planned_arms=30,completed_arms=sum(a['status']=='completed' for a in arms),failed_arms=sum(a['status']=='failed' for a in arms),reference_qualified=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,choices=[0,1],required=True);a=p.parse_args();execute(a.index)
