"""Acquisition-axis photon pooling; post-hoc exposed-data repair, not validation."""
import argparse
import subprocess
import time
import traceback
import numpy as np
from skimage.registration import phase_cross_correlation
import recovery as R
from condition import HERE,ROOT,read,write,sha

DEST=HERE/'structured-repair'

def fit(images,support,rounds):
    assert images.shape==(4*rounds,128,128) and support.shape==images.shape
    summed=images.reshape(rounds,4,128,128).sum(axis=1)
    a=np.array([phase_cross_correlation(summed[0,24:104,24:104],im[24:104,24:104],upsample_factor=50)[0] for im in summed])
    a-=a.mean(axis=0);b=np.zeros((4,2));trace=[]
    def update(groups,group_support,current):
        norm=np.sqrt(np.mean(groups[:,24:104,24:104]**2,axis=(1,2)))
        assert np.all(norm>0),'Uninformative aggregate'
        groups=groups/norm[:,None,None]
        assert np.all(R.warp(group_support,current)[:,24:104,24:104]>=1-1e-6),'Unobserved fit pixels'
        template=R.warp(groups,current).mean(axis=0)[24:104,24:104]
        result=np.array([R.bounded(template,im[24:104,24:104]) for im in groups]);result-=result.mean(axis=0)
        return result
    for step in range(10):
        other=np.tile(b,(rounds,1));aligned=R.warp(images,other).reshape(rounds,4,128,128)
        valid=R.warp(support,other).reshape(rounds,4,128,128)
        a=update(aligned.sum(axis=1),valid.min(axis=1),a)
        other=np.repeat(a,4,axis=0);aligned=R.warp(images,other).reshape(rounds,4,128,128)
        valid=R.warp(support,other).reshape(rounds,4,128,128)
        b=update(aligned.sum(axis=0),valid.min(axis=0),b)
        trace.append(dict(step=step,round_maps=a.tolist(),color_maps=b.tolist()))
    maps=np.repeat(a,4,axis=0)+np.tile(b,(rounds,1))
    assert np.max(np.abs(maps))<24 and np.allclose(maps.mean(axis=0),0)
    assert np.all(R.warp(support,maps)[:,24:104,24:104]>=1-1e-6)
    return maps,trace

def controls():
    from scipy.ndimage import gaussian_filter
    rng=np.random.default_rng(606);records=[]
    for rounds in [7,10]:
        images=np.zeros((rounds,4,128,128))
        for gene in range(60):
            word=np.tile(np.arange(4),3)[:rounds].copy();rng.shuffle(word)
            y,x=rng.integers(32,96,size=2);amplitude=rng.uniform(1,3)
            for r,c in enumerate(word):images[r,c,y,x]+=amplitude
        images=gaussian_filter(images,sigma=(0,0,1.2,1.2)).reshape(4*rounds,128,128)
        rs=np.r_[np.tile([4,2],(rounds//2,1)),np.tile([-4,-2],(rounds//2,1)),np.zeros((rounds%2,2))]
        cs=np.array([[3,2],[1,-1],[-3,-2],[-1,1]])
        for kind in ['zero','round','color','combined']:
            imposed=np.zeros((4*rounds,2))
            if kind in ['round','combined']:imposed+=np.repeat(rs,4,axis=0)
            if kind in ['color','combined']:imposed+=np.tile(cs,(rounds,1))
            observed=R.warp(images,imposed);support=R.warp(np.ones_like(images),imposed)
            maps,trace=fit(observed,support,rounds);applied=np.rint(maps).astype(int)
            exact=np.array_equal(R.warp(observed,applied)[:,16:112,16:112],images[:,16:112,16:112])
            records.append(dict(rounds=rounds,kind=kind,exact_native_image_and_identity_preservation=exact,
                residual_added_rms=float(np.sqrt(np.mean(np.sum((imposed+applied)**2,axis=1)))),maps=maps.tolist()))
    out=dict(passed=all(r['exact_native_image_and_identity_preservation'] for r in records),records=records,seed=606,script_sha256=sha(__file__))
    write(DEST/'controls.json',out);print('CONTROLS',out['passed'],[(r['rounds'],r['kind'],r['residual_added_rms']) for r in records],flush=True)
    assert out['passed'],'Synthetic controls fail; no public repair scoring'

def execute(index):
    for name in ['STRUCTURED_REPAIR.md','structured_repair.py']:
        assert subprocess.check_output(['git','show','HEAD:outputs/barseq3-native-repair-oct05/'+name],cwd=ROOT)==(HERE/name).read_bytes()
    assert read(DEST/'controls.json')['passed'] and read(DEST/'controls.json')['script_sha256']==sha(__file__)
    diagnostic=read(HERE/'objective-diagnosis/field-1/color.json')
    assert diagnostic['objectives']['oracle']['padded']<diagnostic['objectives']['full_model']['padded']
    assert all(d['metrics']['5']['gene_class_identity_f1']>=.9 for d in diagnostic['decoded'].values()),'Initializer repair go/no-go failed'
    book,blank,labels=R.qualify.codebook();names=list(dict.fromkeys(labels));gene_ids=np.array([names.index(n) for n in labels])
    raw,support=R.inputs(index);dest=DEST/f'field-{index}';dest.mkdir(parents=True,exist_ok=True)
    primary=HERE/'recovery'/f'field-{index}'
    for kind in ['zero','common','color']:
        observed=R.warp(raw,R.offsets(kind));observed_support=R.warp(support,R.offsets(kind))
        gp=dest/f'{kind}-geometry.json'
        if gp.exists():geo=read(gp);assert geo['script_sha256']==sha(__file__)
        else:
            start=time.monotonic()
            try:
                maps,trace=fit(observed,observed_support,10);structured_s=time.monotonic()-start
                model=R.qualify.bardensr.registration.find_translations_using_model(observed[:,None].astype(float),book,
                    maximum_wiggle=10,niter=50,initial_guess=np.c_[np.zeros(40),-maps])
                modelmaps=-np.asarray(model.corrections)[:,1:]
                geo=dict(status='completed',maps=dict(structured=maps.tolist(),structured_author=modelmaps.tolist()),trace=trace,
                    author_losses=np.asarray(model.losses).tolist(),structured_s=structured_s,author_s=time.monotonic()-start-structured_s)
            except Exception:geo=dict(status='failed',traceback=traceback.format_exc())
            geo.update(script_sha256=sha(__file__),specification_sha256=sha(HERE/'STRUCTURED_REPAIR.md'));write(gp,geo)
        for policy in ['structured','structured_author']:
            for decoder in R.DECODERS:
                path=dest/f'{kind}__{policy}__{decoder}.json'
                if path.exists():continue
                try:
                    assert geo['status']=='completed'
                    applied=np.rint(geo['maps'][policy]).astype(int);aligned=R.warp(observed,applied)
                    fraction=float(np.all(R.warp(observed_support,applied)[:,24:104,24:104]>=1-1e-6,axis=0).mean());assert fraction>=.99
                    reuse=None
                    for previous in [primary,HERE/'fit-support'/f'field-{index}']:
                        for other in R.POLICIES:
                            source=previous/f'{kind}__{other}__{decoder}.json';r=read(source)
                            if r['status']!='completed':continue
                            if np.array_equal(R.warp(observed,np.array(r['applied_integer_maps']))[:,16:112,16:112],aligned[:,16:112,16:112]):
                                xy=np.array(r['xy']).reshape(-1,2);words=np.array(r['words'],int);runtime=r['runtime']
                                reuse=dict(path=str(source.relative_to(ROOT)),sha256=sha(source),exact_input_equality=True);break
                        if reuse:break
                    if reuse is None:xy,words,runtime=R.decode(aligned,book,decoder)
                    ref=read(primary/f'reference-{decoder}.json')
                    scores={str(r):R.metrics(np.array(ref['xy']).reshape(-1,2),np.array(ref['words'],int),xy,words,gene_ids,names,r) for r in [1.5,3,5]}
                    out=dict(status='completed',source_index=index,kind=kind,policy=policy,decoder=decoder,xy=xy.tolist(),words=words.tolist(),
                        metrics=scores,quality=R.quality(words,blank,fraction),runtime=runtime,reuse=reuse,applied_integer_maps=applied.tolist(),
                        residual_added_rms=float(np.sqrt(np.mean(np.sum((R.offsets(kind)+applied)**2,axis=1)))),geometry_sha256=sha(gp))
                except Exception:out=dict(status='failed',source_index=index,kind=kind,policy=policy,decoder=decoder,traceback=traceback.format_exc())
                out.update(script_sha256=sha(__file__),specification_sha256=sha(HERE/'STRUCTURED_REPAIR.md'));write(path,out)
                print(index,kind,policy,decoder,out['status'],out.get('residual_added_rms'),out.get('metrics',{}).get('5',{}).get('gene_class_identity_f1'),flush=True)
    records=[read(p) for p in DEST.glob('field-*/*__*.json')]
    write(DEST/'state.json',dict(stage='complete' if len(records)==24 else 'partial',planned=24,completed=sum(r['status']=='completed' for r in records),failed=sum(r['status']=='failed' for r in records)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--controls',action='store_true');p.add_argument('--index',type=int,choices=[0,1]);args=p.parse_args()
    if args.controls:controls()
    else:assert args.index is not None;execute(args.index)
