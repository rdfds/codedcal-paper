"""Finite observable-only global search of the pinned single-code objective."""
import argparse
import gc
import subprocess
import time
import traceback
import numpy as np
from scipy.ndimage import gaussian_filter
import recovery as R
from condition import HERE,ROOT,read,write,sha
from bardensr.registration import lowrankregistration as LR

DEST=HERE/'global-repair'
PARTITIONS=[(0,1),(0,2),(0,3)]

def search(images,book,common):
    rounds=images.shape[0]//4
    assert images.shape==(4*rounds,128,128)
    normalized=(book/np.sqrt((book**2).sum(axis=0,keepdims=True))).astype(np.float32)
    denominator=float(np.sum(images.astype(float)**2));records=[]
    assert np.max(np.abs(common))<=10
    # Preserve photons while composing maps. Cropping to 128 between corrections
    # would make the two initialization families have different denominators.
    shifted=R.warp(np.pad(images,((0,0),(30,30),(30,30))),common)
    buffer=np.empty((book.shape[1],168,168),np.float32)
    bestimage=np.empty((168,168),np.float32)
    best=None
    for partition_id,partition in enumerate(PARTITIONS):
        mask=np.isin(np.arange(4*rounds)%4,partition)
        plus=(normalized[mask].T@shifted[mask].reshape(mask.sum(),-1).astype(np.float32)).reshape(-1,188,188)
        minus=(normalized[~mask].T@shifted[~mask].reshape((~mask).sum(),-1).astype(np.float32)).reshape(-1,188,188)
        def evaluate(y,x):
            np.add(plus[:,10-y:178-y,10-x:178-x],minus[:,10+y:178+y,10+x:178+x],out=buffer)
            np.max(buffer,axis=0,out=bestimage)
            loss=1-float(np.sum(bestimage.astype(float)**2))/denominator
            record=dict(partition=partition_id,plus_colors=list(partition),shift_yx=[int(y),int(x)],loss=loss)
            records.append(record)
            return (loss,y*y+x*x,partition_id,int(y),int(x))
        local=min(evaluate(y,x) for y in range(-8,9,2) for x in range(-8,9,2))
        # Each partition gets the same fixed local refinement budget.
        y,x=local[-2:];local=min([local]+[evaluate(y+dy,x+dx) for dy in [-1,0,1] for dx in [-1,0,1]])
        if best is None or local<best:best=local
        del plus,minus;gc.collect()
    partition=PARTITIONS[best[2]];sign=np.where(np.isin(np.arange(4*rounds)%4,partition),1,-1)
    maps=common+sign[:,None]*np.array(best[-2:])[None,:]
    actual=float(LR.calc_loss(images.astype(float),book,-20-maps,zero_padding=20))
    assert abs(actual-best[0])<2e-6,(actual,best[0])
    return maps,dict(grid=records,selected_partition=best[2],selected_plus_colors=list(partition),selected_shift_yx=list(best[-2:]),
                     projected_loss=best[0],actual_author_loss=actual,numerical_loss_error=abs(actual-best[0]))

def common_fit(images,rounds):
    summed=images.reshape(rounds,4,128,128).sum(axis=1)
    initial=np.array([R.bounded(summed[0],im) for im in summed]);initial-=initial.mean(axis=0)
    return np.repeat(np.rint(R.refine(summed,initial)).astype(int),4,axis=0)

def fit(images,support,book):
    start=time.monotonic();common=common_fit(images,images.shape[0]//4)
    families=[]
    for name,initial in [('zero',np.zeros_like(common)),('round',common)]:
        maps,diagnostics=search(images,book,initial)
        families.append(dict(name=name,maps=maps,diagnostics=diagnostics))
    selected=min(families,key=lambda f:(f['diagnostics']['projected_loss'],float(np.sum(f['maps']**2)),f['name']))
    maps=selected['maps'];diagnostics=selected['diagnostics']
    diagnostics.update(selected_round_initialization=selected['name'],families=[dict(name=f['name'],maps=f['maps'].tolist(),diagnostics=f['diagnostics'].copy()) for f in families])
    assert np.all(R.warp(support,maps)[:,24:104,24:104]>=1-1e-6)
    return maps,diagnostics|dict(common_maps=common.tolist(),runtime_s=time.monotonic()-start)

def controls():
    rng=np.random.default_rng(607);records=[]
    for rounds in [7,10]:
        images=np.zeros((rounds,4,128,128));book=np.zeros((rounds,4,60))
        for j in range(60):
            word=np.tile(np.arange(4),3)[:rounds].copy();rng.shuffle(word)
            y,x=rng.integers(32,96,size=2);amplitude=rng.uniform(1,3)
            for r,c in enumerate(word):images[r,c,y,x]+=amplitude;book[r,c,j]=1
        images=gaussian_filter(images,sigma=(0,0,1.2,1.2)).reshape(rounds*4,128,128);book=book.reshape(rounds*4,-1)
        cases=[('zero',np.zeros((rounds*4,2)))]
        rs=np.r_[np.tile([4,2],(rounds//2,1)),np.tile([-4,-2],(rounds//2,1)),np.zeros((rounds%2,2))]
        cases.append(('round',np.repeat(rs,4,axis=0)))
        for p,partition in enumerate(PARTITIONS):
            sign=np.where(np.isin(np.arange(rounds*4)%4,partition),1,-1)
            for amplitude in [(3,1),(4,2),(5,3)]:cases.append((f'partition-{p}-shift-{amplitude[0]}-{amplitude[1]}',sign[:,None]*amplitude))
        for kind,imposed in cases:
            try:
                observed=R.warp(images,imposed);maps,diagnostics=fit(observed,R.warp(np.ones_like(images),imposed),book)
                exact=np.array_equal(R.warp(observed,maps)[:,16:112,16:112],images[:,16:112,16:112])
                record=dict(status='completed',rounds=rounds,kind=kind,exact=exact,residual_added_rms=float(np.sqrt(np.mean(np.sum((imposed+maps)**2,axis=1)))),diagnostics=diagnostics)
            except Exception:record=dict(status='failed',rounds=rounds,kind=kind,exact=False,traceback=traceback.format_exc())
            records.append(record);write(DEST/'controls.json',dict(passed=len(records)==22 and all(r['exact'] for r in records),planned=22,records=records,seed=607,script_sha256=sha(__file__)))
            print('CONTROL',rounds,kind,record['status'],record.get('exact'),record.get('residual_added_rms'),flush=True)
    assert all(r['exact'] for r in records),'Global search control failed; no public scoring'

def execute(index):
    for name in ['GLOBAL_REPAIR.md','global_repair.py']:
        assert subprocess.check_output(['git','show','HEAD:outputs/barseq3-native-repair-oct05/'+name],cwd=ROOT)==(HERE/name).read_bytes()
    assert read(DEST/'controls.json')['passed'] and read(DEST/'controls.json')['script_sha256']==sha(__file__)
    book,blank,labels=R.qualify.codebook();names=list(dict.fromkeys(labels));gene_ids=np.array([names.index(n) for n in labels])
    raw,support=R.inputs(index);dest=DEST/f'field-{index}';dest.mkdir(parents=True,exist_ok=True);primary=HERE/'recovery'/f'field-{index}'
    for kind in ['zero','common','color']:
        imposed=R.offsets(kind);observed=R.warp(raw,imposed);observed_support=R.warp(support,imposed);gp=dest/f'{kind}-geometry.json'
        if gp.exists():geo=read(gp);assert geo['script_sha256']==sha(__file__)
        else:
            try:
                maps,diagnostics=fit(observed,observed_support,book)
                model=R.qualify.bardensr.registration.find_translations_using_model(observed[:,None].astype(float),book,maximum_wiggle=10,niter=50,initial_guess=np.c_[np.zeros(40),-maps])
                geo=dict(status='completed',maps=dict(global_group=maps.tolist(),global_author=(-np.asarray(model.corrections)[:,1:]).tolist()),search=diagnostics,author_losses=np.asarray(model.losses).tolist())
                print('SEARCH',index,kind,diagnostics['selected_partition'],diagnostics['selected_shift_yx'],flush=True)
            except Exception:geo=dict(status='failed',traceback=traceback.format_exc())
            geo.update(script_sha256=sha(__file__),specification_sha256=sha(HERE/'GLOBAL_REPAIR.md'));write(gp,geo)
        for policy in ['global_group','global_author']:
            for decoder in R.DECODERS:
                path=dest/f'{kind}__{policy}__{decoder}.json'
                if path.exists():continue
                try:
                    assert geo['status']=='completed';applied=np.rint(geo['maps'][policy]).astype(int);aligned=R.warp(observed,applied)
                    fraction=float(np.all(R.warp(observed_support,applied)[:,24:104,24:104]>=1-1e-6,axis=0).mean());assert fraction>=.99
                    reuse=None
                    for cohort in ['recovery','fit-support','structured-repair']:
                        policies=R.POLICIES if cohort!='structured-repair' else ['structured','structured_author']
                        for other in policies:
                            source=HERE/cohort/f'field-{index}'/f'{kind}__{other}__{decoder}.json';r=read(source)
                            if r['status']!='completed':continue
                            if np.array_equal(R.warp(observed,np.array(r['applied_integer_maps']))[:,16:112,16:112],aligned[:,16:112,16:112]):
                                xy=np.array(r['xy']).reshape(-1,2);words=np.array(r['words'],int);runtime=r['runtime'];reuse=dict(path=str(source.relative_to(ROOT)),sha256=sha(source),exact_input_equality=True);break
                        if reuse:break
                    if reuse is None:xy,words,runtime=R.decode(aligned,book,decoder)
                    ref=read(primary/f'reference-{decoder}.json');scores={str(r):R.metrics(np.array(ref['xy']).reshape(-1,2),np.array(ref['words'],int),xy,words,gene_ids,names,r) for r in [1.5,3,5]}
                    out=dict(status='completed',source_index=index,kind=kind,policy=policy,decoder=decoder,xy=xy.tolist(),words=words.tolist(),metrics=scores,
                        quality=R.quality(words,blank,fraction),runtime=runtime,reuse=reuse,applied_integer_maps=applied.tolist(),
                        residual_added_rms=float(np.sqrt(np.mean(np.sum((imposed+applied)**2,axis=1)))),geometry_sha256=sha(gp))
                except Exception:out=dict(status='failed',source_index=index,kind=kind,policy=policy,decoder=decoder,traceback=traceback.format_exc())
                out.update(script_sha256=sha(__file__),specification_sha256=sha(HERE/'GLOBAL_REPAIR.md'));write(path,out)
                print(index,kind,policy,decoder,out['status'],out.get('residual_added_rms'),out.get('metrics',{}).get('5',{}).get('gene_class_identity_f1'),flush=True)
    records=[read(p) for p in DEST.glob('field-*/*__*.json')]
    write(DEST/'state.json',dict(stage='complete' if len(records)==24 else 'partial',planned=24,completed=sum(r['status']=='completed' for r in records),failed=sum(r['status']=='failed' for r in records)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--controls',action='store_true');p.add_argument('--index',type=int,choices=[0,1]);args=p.parse_args()
    if args.controls:controls()
    else:assert args.index is not None;execute(args.index)
