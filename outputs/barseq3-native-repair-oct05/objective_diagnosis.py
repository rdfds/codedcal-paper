"""Post-hoc oracle-start controls, without modifying frozen comparison scripts."""
import time
import traceback
import subprocess
import numpy as np
import recovery as R
from condition import HERE,ROOT,read,write,sha
from bardensr.registration import lowrankregistration as LR
from bardensr.registration import translations_tf as T
from bardensr.registration import lowrankregistration_tf as LTF
import tensorflow as tf

DEST=HERE/'objective-diagnosis'

def objective(images,book,maps):
    # An applied image shift m corresponds to source slice origin -m.
    images=np.asarray(images,dtype=float);maps=np.asarray(maps,dtype=float)
    padded=float(LR.calc_loss(images,book,-10-maps,zero_padding=10))
    codes=book/np.sqrt((book**2).sum(axis=0,keepdims=True))
    translated=T.floating_slices(tf.constant(images),tf.constant(24-maps),tf.constant([80,80]),'hermite').numpy()
    dots=np.einsum('fj,fyx->yxj',codes,translated,optimize=True)
    fixed_den=float(np.sum(images[:,24:104,24:104]**2))
    explained=float(np.sum(np.max(dots,axis=-1)**2))
    common=1-explained/fixed_den
    # Also report residual relative to energy at this geometry, to reveal any
    # geometric interpolation-energy confound in fixed-denominator comparison.
    translated_den=float(np.sum(translated**2))
    residual=1-explained/translated_den
    if np.max(np.abs(maps-np.rint(maps)))<1e-10:
        shifted=R.warp(images,maps)[:,24:104,24:104]
        assert np.allclose(translated,shifted,rtol=1e-10,atol=1e-12)
    return dict(padded=padded,fixed_support_fixed_denominator=common,
                fixed_support_geometry_denominator=residual,
                translated_energy_ratio=translated_den/fixed_den)

def execute(index):
    for name in ['OBJECTIVE_DIAGNOSIS.md','objective_diagnosis.py']:
        assert subprocess.check_output(['git','show','HEAD:outputs/barseq3-native-repair-oct05/'+name],cwd=ROOT)==(HERE/name).read_bytes()
    book,blank,labels=R.qualify.codebook();names=list(dict.fromkeys(labels));gene_ids=np.array([names.index(n) for n in labels])
    raw,support=R.inputs(index);dest=DEST/f'field-{index}';dest.mkdir(parents=True,exist_ok=True)
    for kind in ['zero','common','color']:
        outpath=dest/f'{kind}.json'
        if outpath.exists():
            assert read(outpath)['script_sha256']==sha(__file__);continue
        start=time.monotonic()
        try:
            imposed=R.offsets(kind);observed=R.warp(raw,imposed)
            assert np.array_equal(R.warp(observed,-imposed)[:,16:112,16:112],raw[:,16:112,16:112])
            gp=HERE/'recovery'/f'field-{index}'/f'{kind}-geometry.json'
            ip=HERE/'fit-support'/f'field-{index}'/f'{kind}-geometry.json'
            maps=dict(none=np.zeros((40,2)),full_model=np.array(read(gp)['maps']['barcode_model']),
                      interior_model=np.array(read(ip)['maps']['barcode_model']),oracle=-imposed)
            scores={k:objective(observed,book,m) for k,m in maps.items()}
            loss,grad=LR.calc_loss_and_grad(observed.astype(float),book,-10+imposed,zero_padding=10)
            print(index,kind,'objective',scores,'oracle_gradient_max',float(np.abs(grad).max()),flush=True)
            result=R.qualify.bardensr.registration.find_translations_using_model(observed[:,None].astype(float),book,
                maximum_wiggle=10,niter=50,initial_guess=np.c_[np.zeros(40),imposed].astype(float))
            oracle_started=-np.asarray(result.corrections)[:,1:]
            assert np.isfinite(oracle_started).all()
            scores['oracle_started']=objective(observed,book,oracle_started)
            profile=[]
            if kind=='color':
                for fraction in np.linspace(0,1,9):
                    candidate=(1-fraction)*maps['full_model']+fraction*maps['oracle']
                    profile.append(dict(fraction_toward_oracle=float(fraction),objective=objective(observed,book,candidate)))
            correction=np.rint(oracle_started).astype(int);aligned=R.warp(observed,correction)
            actualsupport=R.warp(R.warp(support,imposed),correction)
            fraction=float(np.all(actualsupport[:,24:104,24:104]>=1-1e-6,axis=0).mean());assert fraction>=.99
            decoded={}
            for decoder in R.DECODERS:
                primary=HERE/'recovery'/f'field-{index}';ref=read(primary/f'reference-{decoder}.json');reuse=None
                for policy in R.POLICIES:
                    p=primary/f'{kind}__{policy}__{decoder}.json';prior=read(p)
                    if prior['status']!='completed':continue
                    other=R.warp(observed,np.array(prior['applied_integer_maps']))
                    if np.array_equal(other[:,16:112,16:112],aligned[:,16:112,16:112]):
                        xy=np.array(prior['xy']).reshape(-1,2);words=np.array(prior['words'],int);runtime=prior['runtime']
                        reuse=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),exact_input_equality=True);break
                if reuse is None:xy,words,runtime=R.decode(aligned,book,decoder)
                metrics={str(r):R.metrics(np.array(ref['xy']).reshape(-1,2),np.array(ref['words'],int),xy,words,gene_ids,names,r) for r in [1.5,3,5]}
                decoded[decoder]=dict(xy=xy.tolist(),words=words.tolist(),metrics=metrics,quality=R.quality(words,blank,fraction),runtime=runtime,reuse=reuse)
            out=dict(status='completed',source_index=index,kind=kind,objectives=scores,oracle_gradient=grad.tolist(),
                oracle_gradient_loss=float(loss),oracle_started_maps=oracle_started.tolist(),oracle_started_losses=np.asarray(result.losses).tolist(),
                applied_integer_maps=correction.tolist(),residual_added_rms=float(np.sqrt(np.mean(np.sum((imposed+correction)**2,axis=1)))),
                fractional_added_rms=float(np.sqrt(np.mean(np.sum((imposed+oracle_started)**2,axis=1)))),profile=profile,decoded=decoded,
                primary_geometry_sha256=sha(gp),interior_geometry_sha256=sha(ip),runtime_s=time.monotonic()-start)
        except Exception:out=dict(status='failed',source_index=index,kind=kind,traceback=traceback.format_exc())
        out.update(script_sha256=sha(__file__),specification_sha256=sha(HERE/'OBJECTIVE_DIAGNOSIS.md'))
        write(outpath,out);print(index,kind,out['status'],out.get('residual_added_rms'),flush=True)
    records=[read(p) for p in DEST.glob('field-*/*.json')]
    write(DEST/'state.json',dict(stage='complete' if len(records)==6 else 'partial',planned=6,completed=sum(r['status']=='completed' for r in records),failed=sum(r['status']=='failed' for r in records)))

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,choices=[0,1],required=True);execute(p.parse_args().index)
