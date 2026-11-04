"""Full-cohort post-hoc check of artificial shift-padding influence."""
import argparse
import time
import traceback
import numpy as np
from skimage.registration import phase_cross_correlation
import recovery as R
from condition import HERE,ROOT,read,write,sha

def geometry(observed,support,book):
    start=time.monotonic();summed=observed.reshape(10,4,128,128).sum(axis=1)
    initial=np.array([phase_cross_correlation(summed[0,24:104,24:104],im[24:104,24:104],upsample_factor=50)[0] for im in summed]);initial-=initial.mean(axis=0)
    def refine(images,guess,supports):
        current=guess.copy()
        for _ in range(10):
            valid=R.warp(supports,current)[:,24:104,24:104]
            assert np.all(valid>=1-1e-6),'Fit template uses unobserved pixels'
            template=R.warp(images,current).mean(axis=0)[24:104,24:104]
            current=np.array([R.bounded(template,im[24:104,24:104]) for im in images]);current-=current.mean(axis=0)
        return current
    common=np.repeat(refine(summed,initial,support.reshape(10,4,128,128).min(axis=1)),4,axis=0)
    frame=refine(observed,np.repeat(initial,4,axis=0),support)
    maps=dict(none=np.zeros((40,2)),common_iterative=common,population_template=frame)
    models={};losses={};errors={}
    assert np.all(support[:,24:104,24:104]>=1-1e-6)
    for name,guess in [('common_round',np.repeat(initial,4,axis=0)),('population_template',frame)]:
        try:
            fit=R.qualify.bardensr.registration.find_translations_using_model(observed[:,None,24:104,24:104].astype(float),book,
                maximum_wiggle=10,niter=50,initial_guess=np.c_[np.zeros(40),-guess].astype(float))
            models[name]=-np.asarray(fit.corrections)[:,1:];losses[name]=np.asarray(fit.losses).tolist()
            assert np.isfinite(models[name]).all() and np.isfinite(losses[name]).all()
        except Exception:errors[name]=traceback.format_exc()
    selected=None
    if models:
        selected=min(models,key=lambda k:losses[k][-1]);maps['barcode_model']=models[selected]
    feature=float(np.sqrt(np.mean(np.sum((frame-common)**2,axis=1))))
    return maps,dict(feature=feature,recommend='population_template' if feature>=2 else 'common_iterative',
        model_losses=losses,model_selected_initial=selected,model_failures=errors,fit_s=time.monotonic()-start)

def execute(index):
    assert read(HERE/'recovery/summary.json')['full_cohort_comparison_supported']
    import subprocess
    for name in ['fit_support.py','FIT_SUPPORT_SENSITIVITY.md']:
        assert subprocess.check_output(['git','show','HEAD:outputs/barseq3-native-repair-oct05/'+name],cwd=ROOT)==(HERE/name).read_bytes()
    book,blank,labels=R.qualify.codebook();names=list(dict.fromkeys(labels));gene_ids=np.array([names.index(n) for n in labels])
    raw,support=R.inputs(index);dest=HERE/'fit-support'/f'field-{index}';dest.mkdir(parents=True,exist_ok=True)
    primary=HERE/'recovery'/f'field-{index}'
    references={d:read(primary/f'reference-{d}.json') for d in R.DECODERS}
    for kind in ['zero','common','color']:
        imposed=R.offsets(kind);observed=R.warp(raw,imposed);observed_support=R.warp(support,imposed)
        gp=dest/f'{kind}-geometry.json'
        if gp.exists():
            geo=read(gp);assert geo['fit_support_sha256']==sha(__file__);maps={p:np.array(m) for p,m in geo['maps'].items()}
        else:
            try:
                maps,info=geometry(observed,observed_support,book);maps['oracle']=-imposed
                geo=info|dict(status='completed',maps={p:m.tolist() for p,m in maps.items()},imposed=imposed.tolist(),fit_support_sha256=sha(__file__))
            except Exception:
                maps=dict(none=np.zeros((40,2)),oracle=-imposed)
                geo=dict(status='failed',maps={p:m.tolist() for p,m in maps.items()},traceback=traceback.format_exc(),fit_support_sha256=sha(__file__))
            write(gp,geo)
        for policy in R.POLICIES:
            for decoder in R.DECODERS:
                path=dest/f'{kind}__{policy}__{decoder}.json'
                if path.exists():continue
                try:
                    assert policy in maps,'Calibration did not produce a valid map'
                    correction=np.rint(maps[policy]).astype(int);aligned=R.warp(observed,correction)
                    supportout=R.warp(observed_support,correction);fraction=float(np.all(supportout[:,24:104,24:104]>=1-1e-6,axis=0).mean())
                    assert fraction>=.99,'Applied correction loses fixed scoring support'
                    reused=None
                    for other in R.POLICIES:
                        candidate=primary/f'{kind}__{other}__{decoder}.json';record=read(candidate)
                        if record['status']!='completed':continue
                        prior=R.warp(observed,np.array(record['applied_integer_maps']))
                        if np.array_equal(prior[:,16:112,16:112],aligned[:,16:112,16:112]):
                            xy=np.array(record['xy'],float).reshape(-1,2);words=np.array(record['words'],int);runtime=record['runtime']|dict(actual_execution_kernel_s=0.)
                            reused=dict(path=str(candidate.relative_to(ROOT)),sha256=sha(candidate),exact_decoder_input_equality=True);break
                    if reused is None:xy,words,runtime=R.decode(aligned,book,decoder)
                    ref=references[decoder]
                    if policy=='oracle':assert xy.tolist()==ref['xy'] and words.tolist()==ref['words']
                    scores={str(r):R.metrics(np.array(ref['xy']).reshape(-1,2),np.array(ref['words'],int),xy,words,gene_ids,names,r) for r in [1.5,3,5]}
                    out=dict(status='completed',kind=kind,policy=policy,decoder=decoder,source_index=index,
                        xy=xy.tolist(),words=words.tolist(),metrics=scores,quality=R.quality(words,blank,fraction),runtime=runtime,
                        reused_primary=reused,applied_integer_maps=correction.tolist(),geometry_sha256=sha(gp),
                        fit_support_sha256=sha(__file__),recovery_sha256=sha(HERE/'recovery.py'),specification_sha256=sha(HERE/'FIT_SUPPORT_SENSITIVITY.md'))
                except Exception:out=dict(status='failed',kind=kind,policy=policy,decoder=decoder,source_index=index,traceback=traceback.format_exc(),fit_support_sha256=sha(__file__))
                write(path,out);print(index,kind,policy,decoder,out['status'],out.get('metrics',{}).get('5',{}).get('gene_class_identity_f1'),flush=True)
    arms=[read(p) for p in dest.glob('*__*.json')];assert len(arms)==30
    write(dest/'completion.json',dict(planned=30,completed=sum(a['status']=='completed' for a in arms),failed=sum(a['status']=='failed' for a in arms)))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,choices=[0,1],required=True);a=p.parse_args();execute(a.index)
