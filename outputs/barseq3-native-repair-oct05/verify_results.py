"""Recount retained calls and validate lineage without dense decoder caches."""
import ast
import sys
import h5py
import numpy as np
from condition import HERE,ROOT,CACHE,read,write,sha
sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal.evaluation import spatial_match

def frozen_function(path,name,namespace):
    node=next(n for n in ast.parse(path.read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),namespace)
    return namespace[name]

def equivalent(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for k in a:equivalent(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):equivalent(x,y)
    elif a is None or isinstance(a,(str,bool)):assert a==b
    else:assert abs(a-b)<1e-12,(a,b)

def main():
    study=HERE.parent/'independent-comparison-oct05'
    score=frozen_function(study/'run.py','score',dict(np=np,spatial_match=spatial_match,inside=lambda xy:np.all((xy>=24)&(xy<176),axis=1)))
    integer_verified=0;recovery_verified=0;files=[]
    for split in ['development','validation']:
        source=read(study/split/'source.json')
        path=ROOT.parent/'bardensr-upstream.nosync/bardensr/benchmarks'/(source['name']+'.hdf5')
        assert sha(path)==source['sha256']
        with h5py.File(path) as f:
            xy=np.c_[f['rolonies/m2'][:],f['rolonies/m1'][:]];words=f['rolonies/j'][:].astype(int)
            good=(f['rolonies/status'][:].astype('U')=='good')&np.all((xy>=24)&(xy<176),axis=1)
            xy,words=xy[good],words[good];nwords=f['codebook'].shape[-1]
        for kind,dose in [('zero',0),('round',2),('round',4),('color',2),('color',4)]:
            p=study/'integer-control'/split/f'{kind}-a{dose}.json';case=read(p)
            assert case['integer_control_sha256']==sha(study/'integer_control.py')
            assert case['specification_sha256']==sha(study/'INTEGER_CONTROL.md')
            assert len(case['arms'])==12
            for arm in case['arms'].values():
                prediction=np.array(arm['xy'],float).reshape(-1,2);ids=np.array(arm['words'],int)
                for radius in [1.5,3,5]:equivalent(score(xy,words,prediction,ids,nwords,radius),arm['metrics'][str(radius)])
                maps=np.array(arm['correction_yx']);assert np.array_equal(maps,np.rint(maps))
                integer_verified+=1
            files.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    assert integer_verified==120
    complete=read(HERE/'recovery/summary.json');assert complete['full_cohort_comparison_supported']
    classes=read(HERE/'codebook-classes.json');labels=classes['names'];blank=np.array(classes['blank'],bool)
    names=list(dict.fromkeys(labels));ids=np.array([names.index(n) for n in labels])
    metrics=frozen_function(HERE/'recovery.py','metrics',dict(np=np,spatial_match=spatial_match))
    for index in [0,1]:
        qualification=read(HERE/f'field-{index}/qualification.json');assert qualification['qualified']
        for r in range(10):
            marker=read(HERE/f'field-{index}/round-{r}-demons.json')
            assert marker['output_sha256']==sha(CACHE/f'native-repair-source-{index}-round-{r}.npz')
            assert marker['condition_sha256']==sha(HERE/'condition.py')
        folder=HERE/'recovery'/f'field-{index}'
        references={d:read(folder/f'reference-{d}.json') for d in ['correlation','regression']}
        common=np.array(read(folder/'common-geometry.json')['imposed']);color=np.array(read(folder/'color-geometry.json')['imposed'])
        assert sorted(map(tuple,common))==sorted(map(tuple,color))
        for kind in ['zero','common','color']:
            gp=folder/f'{kind}-geometry.json'
            for policy in ['none','common_iterative','population_template','barcode_model','oracle']:
                for decoder in ['correlation','regression']:
                    p=folder/f'{kind}__{policy}__{decoder}.json';arm=read(p);ref=references[decoder]
                    assert arm['status']=='completed' and arm['recovery_sha256']==sha(HERE/'recovery.py')
                    assert arm['geometry_sha256']==sha(gp) and arm['reference_sha256']==sha(folder/f'reference-{decoder}.json')
                    xy=np.array(arm['xy'],float).reshape(-1,2);words=np.array(arm['words'],int)
                    for radius in [1.5,3,5]:
                        equivalent(metrics(np.array(ref['xy']).reshape(-1,2),np.array(ref['words'],int),xy,words,ids,names,radius),arm['metrics'][str(radius)])
                    nblank=int(blank[words].sum());ngene=len(words)-nblank
                    assert arm['quality']['nonblank']==ngene and arm['quality']['blanks']==nblank
                    assert abs(arm['quality']['proxy']-nblank/max(ngene,1)*int((~blank).sum())/int(blank.sum()))<1e-12
                    if policy=='oracle':assert arm['xy']==ref['xy'] and arm['words']==ref['words']
                    maps=np.array(arm['applied_integer_maps']);assert np.array_equal(maps,np.rint(maps))
                    recovery_verified+=1;files.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p)))
    assert recovery_verified==60
    write(HERE/'verification.json',dict(integer_arms_verified=integer_verified,recovery_arms_verified=recovery_verified,
        radii=[1.5,3,5],native_array_hashes_verified=20,matched_full_offset_distributions_verified=True,
        note='Independent retained-call recount and lineage check using frozen scoring definitions; not independent biological validation.',
        scripts_sha256=sha(__file__),files=files))
    print('Recounted 120 integer and 60 recovery arms at every radius; all lineage checks pass',flush=True)

if __name__=='__main__':main()
