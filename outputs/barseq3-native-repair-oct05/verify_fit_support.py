"""Verify every interior-fit arm, including failures and exact input reuse."""
import numpy as np
from scipy.ndimage import shift
from condition import HERE,ROOT,CACHE,read,write,sha
from verify_results import frozen_function,equivalent,spatial_match

def warp(images,maps):
    return np.stack([shift(im,t,order=1,mode='constant',cval=0,prefilter=False) for im,t in zip(images,maps)])

def main():
    summary=read(HERE/'fit-support/summary.json')
    assert summary['counts']['missing']==0,'Incomplete follow-up cannot be verified as complete'
    classes=read(HERE/'codebook-classes.json');labels=classes['names'];blank=np.array(classes['blank'],bool)
    names=list(dict.fromkeys(labels));ids=np.array([names.index(n) for n in labels])
    metrics=frozen_function(HERE/'recovery.py','metrics',dict(np=np,spatial_match=spatial_match))
    completed=failed=reused=0;files=[]
    for index in [0,1]:
        arrays=[];supports=[]
        for r in range(10):
            p=CACHE/f'native-repair-source-{index}-round-{r}.npz'
            assert sha(p)==read(HERE/f'field-{index}/round-{r}-demons.json')['output_sha256']
            with np.load(p) as a:arrays.extend(a['images']);supports.extend(a['support'])
        raw=np.stack(arrays);support=np.stack(supports)
        folder=HERE/'fit-support'/f'field-{index}';primary=HERE/'recovery'/f'field-{index}'
        for kind in ['zero','common','color']:
            gp=folder/f'{kind}-geometry.json';geo=read(gp)
            assert geo['fit_support_sha256']==sha(HERE/'fit_support.py')
            imposed=np.array(read(primary/f'{kind}-geometry.json')['imposed'])
            observed=warp(raw,imposed);observed_support=warp(support,imposed)
            for policy in ['none','common_iterative','population_template','barcode_model','oracle']:
                for decoder in ['correlation','regression']:
                    p=folder/f'{kind}__{policy}__{decoder}.json';arm=read(p)
                    assert arm['fit_support_sha256']==sha(HERE/'fit_support.py')
                    files.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),status=arm['status']))
                    if arm['status']=='failed':
                        assert arm.get('traceback');failed+=1;continue
                    assert arm['status']=='completed'
                    assert arm['recovery_sha256']==sha(HERE/'recovery.py')
                    assert arm['specification_sha256']==sha(HERE/'FIT_SUPPORT_SENSITIVITY.md')
                    assert arm['geometry_sha256']==sha(gp)
                    maps=np.array(arm['applied_integer_maps']);assert np.array_equal(maps,np.rint(maps))
                    assert np.array_equal(maps,np.rint(np.array(geo['maps'][policy])))
                    fraction=float(np.all(warp(observed_support,maps)[:,24:104,24:104]>=1-1e-6,axis=0).mean())
                    assert fraction==arm['quality']['support']
                    ref=read(primary/f'reference-{decoder}.json')
                    xy=np.array(arm['xy'],float).reshape(-1,2);words=np.array(arm['words'],int)
                    for radius in [1.5,3,5]:
                        equivalent(metrics(np.array(ref['xy']).reshape(-1,2),np.array(ref['words'],int),xy,words,ids,names,radius),arm['metrics'][str(radius)])
                    ngene=int((~blank[words]).sum());nblank=len(words)-ngene
                    assert ngene==arm['quality']['nonblank'] and nblank==arm['quality']['blanks']
                    proxy=nblank/max(ngene,1)*int((~blank).sum())/int(blank.sum())
                    assert abs(proxy-arm['quality']['proxy'])<1e-12
                    if policy=='oracle':assert arm['xy']==ref['xy'] and arm['words']==ref['words']
                    reuse=arm['reused_primary']
                    if reuse:
                        source=ROOT/reuse['path'];assert sha(source)==reuse['sha256']
                        prior=read(source)
                        assert arm['xy']==prior['xy'] and arm['words']==prior['words']
                        assert np.array_equal(warp(observed,maps)[:,16:112,16:112],warp(observed,np.array(prior['applied_integer_maps']))[:,16:112,16:112])
                        assert arm['runtime']['actual_execution_kernel_s']==0 and reuse['exact_decoder_input_equality']
                        reused+=1
                    completed+=1
    assert completed+failed==60
    assert completed==summary['counts']['completed'] and failed==summary['counts']['failed']
    write(HERE/'fit-support/verification.json',dict(completed_verified=completed,failures_retained=failed,exact_input_reuses_verified=reused,
        radii=[1.5,3,5],scripts_sha256=sha(__file__),files=files,
        scope='Retained-call recount, lineage and exact reuse verification; does not provide independent biological truth'))
    print('Verified',completed,'completed arms;',failed,'retained failures;',reused,'exact input reuses',flush=True)

if __name__=='__main__':main()
