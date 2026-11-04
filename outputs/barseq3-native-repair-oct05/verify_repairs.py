"""Retained-call, photon-support and whole-domain assignment checks for repairs."""
import numpy as np
from scipy.ndimage import shift
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from condition import HERE,ROOT,CACHE,read,write,sha
from verify_results import frozen_function,equivalent,spatial_match

def main():
    classes=read(HERE/'codebook-classes.json');names=list(dict.fromkeys(classes['names']))
    assert classes['codebook_sha256']==sha(HERE.parent/'independent-comparison-oct05/barseq3-adult-codebook.csv')
    gene_ids=np.array([names.index(n) for n in classes['names']]);blank=np.array(classes['blank'],bool)
    book=np.zeros((40,len(blank)))
    for j,s in enumerate(classes['sequences']):
        for r,c in enumerate(s):book[4*r+'GTCA'.index(c),j]=1
    normalized=book/np.sqrt((book**2).sum(axis=0,keepdims=True))
    metrics=frozen_function(HERE/'recovery.py','metrics',dict(np=np,spatial_match=spatial_match))
    warp=frozen_function(HERE/'recovery.py','warp',dict(np=np,shift=shift))
    offsets=frozen_function(HERE/'recovery.py','offsets',dict(np=np))
    files=[];assignment=[];checks=[];arm_count=0;arrays=0;direct_loss=[]
    for index in [0,1]:
        images=[];supports=[]
        for r in range(10):
            path=CACHE/f'native-repair-source-{index}-round-{r}.npz';marker=read(HERE/f'field-{index}/round-{r}-demons.json')
            assert sha(path)==marker['output_sha256'] and marker['condition_sha256']==sha(HERE/'condition.py')
            with np.load(path) as d:images.extend(d['images']);supports.extend(d['support'])
            arrays+=1
        images=np.stack(images);supports=np.stack(supports)
        for kind in ['zero','common','color']:
            observed=warp(images,offsets(kind));support=warp(supports,offsets(kind))
            for cohort,script,protocol in [('objective-diagnosis','objective_diagnosis.py','OBJECTIVE_DIAGNOSIS.md'),
                                           ('structured-repair','structured_repair.py','STRUCTURED_REPAIR.md'),
                                           ('global-repair','global_repair.py','GLOBAL_REPAIR.md')]:
                folder=HERE/cohort/f'field-{index}'
                if cohort=='objective-diagnosis':
                    path=folder/f'{kind}.json';record=read(path);assert record['status']=='completed'
                    files.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path)))
                    records=[dict(record,decoder=d,policy='oracle_started',**{k:v for k,v in arm.items() if k not in record}) for d,arm in record['decoded'].items()]
                else:
                    paths=list(folder.glob(kind+'__*.json'));assert len(paths)==4
                    records=[]
                    for path in paths:
                        record=read(path);assert record['geometry_sha256']==sha(folder/f'{kind}-geometry.json');records.append(record)
                        files.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path)))
                for arm in records:
                    assert arm['status']=='completed' and arm['script_sha256']==sha(HERE/script)
                    assert arm['specification_sha256']==sha(HERE/protocol)
                    applied=np.array(arm['applied_integer_maps']);assert np.array_equal(applied,np.rint(applied))
                    aligned=warp(observed,applied);fraction=float(np.all(warp(support,applied)[:,24:104,24:104]>=1-1e-6,axis=0).mean())
                    assert abs(fraction-arm['quality']['support'])<1e-12
                    reuse=arm.get('reuse')
                    if reuse:
                        prior=ROOT/reuse['path'];assert sha(prior)==reuse['sha256'];old=read(prior)
                        assert np.array_equal(warp(observed,np.array(old['applied_integer_maps']))[:,16:112,16:112],aligned[:,16:112,16:112])
                        assert arm['xy']==old['xy'] and arm['words']==old['words']
                    ref=read(HERE/'recovery'/f'field-{index}'/f"reference-{arm['decoder']}.json")
                    a=np.array(ref['xy']).reshape(-1,2);b=np.array(arm['xy']).reshape(-1,2);words=np.array(arm['words'],int)
                    for radius in [1.5,3,5]:
                        expected=metrics(a,np.array(ref['words'],int),b,words,gene_ids,names,radius)
                        equivalent(expected,arm['metrics'][str(radius)])
                        dist=cdist(a,b);edges=dist<=radius;n,m=dist.shape
                        cardinality=int((maximum_bipartite_matching(csr_matrix(edges),perm_type='column')>=0).sum()) if n and m else 0
                        pairs,diag=spatial_match(a,b,radius);assert len(pairs)==cardinality and diag['skipped_components']==0
                        reward=(min(n,m)+1)*radius;cost=np.zeros((n,m+n));cost[:,:m]=np.where(edges,dist-reward,reward)
                        rr,cc=linear_sum_assignment(cost);good=cc<m;rr,cc=rr[good],cc[good];good=edges[rr,cc];rr,cc=rr[good],cc[good]
                        assert len(rr)==cardinality
                        assert abs(float(dist[rr,cc].sum())-float(np.linalg.norm(a[pairs[:,0]]-b[pairs[:,1]],axis=1).sum()))<1e-8
                        correct=int((gene_ids[np.array(ref['words'],int)[rr]]==gene_ids[words[cc]]).sum())
                        assignment.append(dict(cohort=cohort,index=index,kind=kind,policy=arm['policy'],decoder=arm['decoder'],radius=radius,
                            primary_f1=expected['gene_class_identity_f1'],global_f1=2*correct/(n+m) if n+m else 0,
                            tied_identity_difference=correct-expected['correct_gene_classes']))
                    nblank=int(blank[words].sum());ngene=len(words)-nblank
                    assert arm['quality']['blanks']==nblank and arm['quality']['nonblank']==ngene
                    assert abs(arm['quality']['proxy']-nblank/max(ngene,1)*int((~blank).sum())/int(blank.sum()))<1e-12
                    residual=float(np.sqrt(np.mean(np.sum((offsets(kind)+applied)**2,axis=1))))
                    assert abs(residual-arm['residual_added_rms'])<1e-12
                    if residual==0:
                        assert np.array_equal(aligned[:,16:112,16:112],images[:,16:112,16:112])
                        assert arm['xy']==ref['xy'] and arm['words']==ref['words']
                    checks.append(dict(cohort=cohort,index=index,kind=kind,policy=arm['policy'],decoder=arm['decoder'],exact_reference=residual==0))
                    arm_count+=1
            # Direct float64 projection of independently shifted full frames;
            # distinct from the accelerated grouped-projection implementation.
            geo=read(HERE/'global-repair'/f'field-{index}'/f'{kind}-geometry.json')
            maps=np.array(geo['maps']['global_group']);full=warp(np.pad(observed,((0,0),(20,20),(20,20))),maps).reshape(40,-1)
            maximum=np.full(full.shape[1],-np.inf)
            for start in range(0,book.shape[1],128):
                dots=normalized[:,start:start+128].T@full;maximum=np.maximum(maximum,dots.max(axis=0))
            actual=1-float(np.sum(maximum**2))/float(np.sum(observed.astype(float)**2));stored=geo['search']['projected_loss']
            assert abs(actual-stored)<2e-6
            direct_loss.append(dict(index=index,kind=kind,direct_numpy_loss=actual,accelerated_loss=stored,error=abs(actual-stored)))
    assert arm_count==60 and arrays==20
    write(HERE/'repair-verification.json',dict(verified_arms=arm_count,radii=[1.5,3,5],native_array_hashes=arrays,checks=checks,
        assignments=assignment,direct_losses=direct_loss,files=files,script_sha256=sha(__file__),
        equal_cost_identity_differences=sum(r['tied_identity_difference']!=0 for r in assignment),
        largest_observed_equal_cost_f1_change=max(abs(r['primary_f1']-r['global_f1']) for r in assignment),
        note='Numerical/provenance verification; not independent biological confirmation. Two equal-cost assignments are not an exhaustive identity bound.'))
    print('Verified 60 repair/control arms, 180 assignments, 20 native arrays and six direct global losses',flush=True)

if __name__=='__main__':main()
