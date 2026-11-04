"""Independent recounts of retained supports, simulation records and hashes."""
import json
import numpy as np
from common import P,ROOT,read,save,sha,cd20_registry,historical_members,input_path
from reference_characterize import registry

def main():
    comparison=read(str((P/'comparison.json').relative_to(ROOT)))
    checks=[];control_count=0;group_count=0;source_count=0
    assert comparison['distinct_acquisitions']==8 and len(comparison['acquisitions'])==9
    old_members=historical_members();primary_fingerprints=[]
    registry_by_tag={r['tag']:r for r in registry()}
    for row in cd20_registry():
        member_set=[]
        for r in range(1,5):
            name=row['prefix']+f'RESI/RESI rounds/RESI_R{r}_apicked.hdf5'
            member_set.append(next(v['sha256'] for k,v in old_members[name].items()
                if k in ('bsdtar','unar','unrar') and 'sha256' in v))
        primary_fingerprints.append(tuple(member_set))
    assert len(set(primary_fingerprints))==6
    for a in comparison['acquisitions']:
        tag=a['tag'];d=np.load(P/f'reference-{tag}.npz')
        assert len(d['that'])==a['passing_sources']
        source_count+=len(d['that'])
        s=d['support'];n=s[:,2:6];loc=s[:,6:10]
        def triple(v):
            return v[:,0]*v[:,1]*(v[:,2]+v[:,3])+v[:,2]*v[:,3]*(v[:,0]+v[:,1])
        assert np.array_equal(triple(n),d['event_countw'])
        assert np.array_equal(triple(loc),d['localization_countw'])
        assert np.array_equal(s[:,0].astype(int),d['round'])
        assert np.array_equal(s[:,1].astype(int),d['pore'])
        for key,cols in [('eligible_per_panel_event_count_quantiles',slice(2,6)),
                         ('eligible_per_panel_localization_count_quantiles',slice(6,10))]:
            assert np.array_equal(np.quantile(s[:,cols],[.1,.5,.9],axis=0),a['empirical_conditions'][key])
        for name,key in [('unweighted',None),('event_count','event_countw'),('localization_count','localization_countw')]:
            estimate=float(np.average(d['that'],weights=d[key] if key else None))
            assert abs(estimate-a['estimates'][name]['raw'])<1e-12
        if a['reference_replay'] is not None:assert a['reference_replay']['passed'],tag
        assert a['full_code_moment_replay']['passed'],(tag,a['full_code_moment_replay'])
        for member in a['input_members']:
            assert sha(input_path(registry_by_tag[tag],member['round']))==member['sha256'],member
            r=member['round']
            if tag=='npc':
                prior=read('outputs/shaw-sensitivity-sep26/input_manifest.json')['raw']['reinhardt']
                name=f'Figure2,EDF06,EDF07/RESI/RESI rounds/RESI_R{r}_apicked.hdf5'
                assert prior[name]['unrar']['sha256']==member['sha256']
            elif tag.startswith('sdc_'):
                prior=read('outputs/heldout-nup96-sep30/sdc-recovery.json')[f'round{r}']
                assert prior['reduced_sha256']==member['sha256']
                assert prior['source_sha256']==member['source_sha256']
            else:
                row=next(r for r in cd20_registry() if r['tag']==tag)
                name=row['prefix']+f'RESI/RESI rounds/RESI_R{r}_apicked.hdf5'
                prior=old_members[name]
                expected=next(v['sha256'] for k,v in prior.items() if k in ('unrar','unar','bsdtar') and 'sha256' in v)
                assert expected==member['sha256']
        if tag in [r['tag'] for r in cd20_registry()]:
            c=read(str((P/f'controls-{tag}.json').relative_to(ROOT)))
            assert c['complete'] and len(c['records'])==16
            assert sha(P/f'empirical-{tag}.npz')==c['metadata']['empirical_sha256']
            observed={(r['model'],r['T'],r['seed']) for r in c['records']}
            assert observed=={(m,t,s) for m in ('gaussian','residual') for t in (0.,1.,3.,5.) for s in (2610071101,2610071102)}
            for r in c['records']:
                group=np.load(P/r['group_records'])['by_group']
                assert group.shape==(4000,6),r
                estimate=group[:,:3].sum(0)/group[:,3:].sum(0)
                assert np.allclose(estimate,r['estimate'],rtol=0,atol=1e-12)
                se=np.sqrt(len(group))*(group[:,:3]-group[:,3:]*estimate).std(0,ddof=1)/group[:,3:].sum(0)
                assert np.allclose(se,r['se'],rtol=0,atol=1e-12)
                assert int(group[:,3].sum())==r['passing_sources']
                group_count+=len(group);control_count+=1
        checks.append(dict(view=tag,source_replay=a['reference_replay'],moments=a['full_code_moment_replay'],recounts_pass=True))
    assert control_count==96 and group_count==384000
    # Entire matched-seed world checks the sparse implementation against a
    # retained dense-implementation primary result, beyond small graph tests.
    from simulate_reference import templates,run_world
    template,qc,_=templates('edf10rep1')
    sparse=run_world(template,qc,0.,'gaussian',2610071101,4000)
    dense=np.load(P/'control-edf10rep1-gaussian-T0-s2610071101.npz')['by_group']
    sparse_difference=float(np.max(abs(sparse['by_group']-dense)))
    assert np.allclose(sparse['by_group'],dense,rtol=1e-10,atol=1e-8),sparse_difference
    hashes={str(f.relative_to(P)):sha(f) for f in sorted(P.rglob('*')) if f.is_file()
            and f.name!='verification.json' and '__pycache__' not in f.parts}
    out=dict(all_pass=True,views=9,distinct_acquisitions=8,unique_cd20_member_sets=6,
             source_rows_recounted=source_count,primary_control_worlds=control_count,
             primary_group_draws=group_count,checks=checks,artifact_sha256=hashes,
             sparse_dense_world=dict(groups=4000,max_abs_difference=sparse_difference,passed=True),
             scope='Computational/provenance verification. Inferred clusters and noise models do not provide physical source truth.')
    save(P/'verification.json',out);print(json.dumps({k:v for k,v in out.items() if k not in ('checks','artifact_sha256')},indent=2))

if __name__=='__main__':main()
