"""Consolidate every acquisition, uncertainty, and conditional reference fit."""
import csv
import importlib.util
import json
import numpy as np
from scipy.stats import t as student
from common import P,ROOT,cd20_registry,read,save,sha
from reference_characterize import registry,load_round

def replay_legacy_code(row):
    # Historical CD20 code phases were evaluated on global float32 coordinates;
    # target-matching diagnostics use centered float64 coordinates. Re-execute
    # the retained historical function, rather than loosening a replay tolerance.
    spec=importlib.util.spec_from_file_location('historical_cd20_code',
        ROOT/'outputs/replication-cd20-sep15/analyze_field.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    centers=np.load(row['center'])['centers'];rounds={}
    for r in range(1,5):
        d=load_round(row,r,centers)
        for key in ('lpx','lpy'):d[key]=d[key].astype(float)
        rounds[r]=d
    return module.code_moment_components(rounds,len(centers),
        {int(k):v for k,v in row['panels'].items()})['pooled']

def replay_legacy_sdc(row):
    # Scalar per-angle operations in the original SDC code retain float32
    # phases even after centering. The vectorized diagnostic promotes them.
    spec=importlib.util.spec_from_file_location('historical_sdc_code',
        ROOT/'outputs/independent-field-sep15/analyze_field.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.SCRATCH=ROOT.parent/'heldout-nup96.nosync/independent-field'
    z,n=module.code_moments(row['tag']=='sdc_corrected')
    return module.split_variants(z,n)['pooled']

def ratio_set(blocks):
    # Blocks are sums of source/group contributions: use equal spatial-block
    # resampling units, rather than treating emitters as independent replicates.
    n=len(blocks);mean=blocks.mean(0);cov=np.cov(blocks.T,ddof=1)/n
    x,y=mean;crit=student.ppf(.975,n-1)**2
    a=y*y-crit*cov[1,1];b=-2*(x*y-crit*cov[0,1]);c=x*x-crit*cov[0,0]
    disc=b*b-4*a*c
    if a<=0 or disc<0:
        return dict(kind='unbounded_or_empty',ratio=float(x/y),denominator_positive=y>0,
                    bounded=False,polynomial=[a,b,c])
    roots=np.sort([(-b-np.sqrt(disc))/(2*a),(-b+np.sqrt(disc))/(2*a)])
    return dict(kind='finite',ratio=float(x/y),interval=roots,bounded=True,
                denominator_positive=y>0,polynomial=[a,b,c])

def spatial_intervals(tag):
    d=np.load(P/f'reference-{tag}.npz');centers=d['centers'][:,:2]
    lo=centers.min(0);extent=np.maximum(centers.max(0)-lo,1.)
    cell=np.minimum(((centers-lo)/extent*[4,3]).astype(int),[3,2]);tile=cell[:,0]*3+cell[:,1]
    out={};block_store={}
    for name in ('full_code','selected_code','diagonal'):
        comp=d[name];blocks=np.stack([comp[tile==t].sum(0) for t in np.unique(tile)])
        block_store[name]=blocks
        per_q=[]
        for j,q in enumerate((.10,.15)):
            x=blocks[:,j,:,0].mean(1);y=blocks[:,j,:,1].mean(1)
            interval=ratio_set(np.column_stack([x,y]));interval['q']=q
            interval['T_at_q']=float(-np.log(interval['ratio'])/q**2) if interval['ratio']>0 and interval['denominator_positive'] else None
            if interval['bounded'] and interval['interval'][0]>0 and interval['denominator_positive']:
                interval['T_interval']=(-np.log(interval['interval'][::-1])/q**2)
            else:interval['T_interval']=None
            per_q.append(interval)
        out[name]=dict(blocks=len(blocks),frequencies=per_q)
    rng=np.random.default_rng(26100761);n=len(np.unique(tile));indices=rng.integers(n,size=(2000,n))
    pair=[]
    for j,q in enumerate((.10,.15)):
        ratios=[];denominators=[]
        for name in ('selected_code','diagonal'):
            b=block_store[name][indices].sum(1)[:,j]
            den=b[:,:,1].mean(1);denominators.append(den)
            ratios.append(b[:,:,0].mean(1)/den)
        valid=(ratios[0]>0)&(ratios[1]>0)&(denominators[0]>0)&(denominators[1]>0)&np.isfinite(ratios[0])&np.isfinite(ratios[1])
        gap=np.log(ratios[1][valid]/ratios[0][valid])/q**2
        estimate=out['selected_code']['frequencies'][j]['T_at_q']
        reference=out['diagonal']['frequencies'][j]['T_at_q']
        pair.append(dict(q=q,code_minus_diagonal=estimate-reference if estimate is not None and reference is not None else None,
            valid_fraction=float(valid.mean()),
            conditional_interval=np.quantile(gap,[.025,.975]) if valid.mean()>=.95 else None,
            note='Paired spatial-block bootstrap; no equivalence margin or physical-truth claim.'))
    out['paired_selected_difference']=pair
    return out

def responses(tag):
    controls=read(str((P/f'controls-{tag}.json').relative_to(ROOT)))
    assert controls['complete'] and len(controls['records'])==16
    result={}
    for index,name in enumerate(('unweighted','event_count','localization_count')):
        fits={}
        for model in ('gaussian','residual'):
            rows=[r for r in controls['records'] if r['model']==model]
            x=np.array([[1.,r['T']] for r in rows]);y=np.array([r['estimate'][index] for r in rows])
            se=np.maximum([r['se'][index] for r in rows],.001);xx=x/se[:,None];yy=y/se
            beta=np.linalg.lstsq(xx,yy,rcond=None)[0];res=y-x@beta
            # The same seed supplies common random draws at all T settings.
            # Preserve that pairing instead of treating eight world means as
            # independent. Different seeds form independent simulation blocks.
            covariance_y=np.zeros((len(rows),len(rows)))
            for seed in (2610071101,2610071102):
                indices=[i for i,r in enumerate(rows) if r['seed']==seed]
                influence=[]
                for i in indices:
                    group=np.load(P/rows[i]['group_records'])['by_group']
                    influence.append((group[:,index]-y[i]*group[:,index+3])/group[:,index+3].mean())
                influence=np.stack(influence,1)
                covariance_y[np.ix_(indices,indices)]=np.cov(influence.T,ddof=1)/len(influence)
            operator=np.linalg.inv(xx.T@xx)@x.T@np.diag(1/se**2)
            covariance=operator@covariance_y@operator.T
            bse=np.sqrt(covariance[1,1]);ci=beta[1]+np.array([-1,1])*student.ppf(.975,3999)*bse
            fits[model]=dict(intercept=beta[0],slope=beta[1],slope_interval=ci,
                             max_affine_residual=float(abs(res).max()),
                             covariance_method='Paired per-group influence across T within each seed; two independent seed blocks. Conditional on inferred templates/noise.',
                             mean_event_coverage=float(np.mean([r['event_coverage'] for r in rows])),
                             mean_impure_fraction=float(np.mean([r['impure_fraction'] for r in rows])))
        slopes=[f['slope'] for f in fits.values()];difference=abs(slopes[0]-slopes[1])/max(abs(np.mean(slopes)),1e-12)
        qualifies=bool(all(f['slope_interval'][0]>0 and f['max_affine_residual']<=.2 for f in fits.values()) and difference<=.1)
        real=read(str((P/f'characterization-{tag}.json').relative_to(ROOT)))['estimates'][name]
        sensitivity={}
        for model,fit in fits.items():
            sensitivity[model]=dict(conditional_T=(real['raw']-fit['intercept'])/fit['slope'],
                transformed_observation_interval=(np.array(real['interval'])-fit['intercept'])/fit['slope'],
                note='Observation interval only; model and response-fit uncertainty are separate, not physical accuracy.')
        result[name]=dict(models=fits,relative_slope_difference=difference,
                          qualified_conditional_response=qualifies,conditional_sensitivity=sensitivity)
    return result

def main():
    out=[];rows=registry();members=__import__('common').historical_members()
    alias=[]
    for r in range(1,5):
        names=[p+f'RESI/RESI rounds/RESI_R{r}_apicked.hdf5' for p in ('Figure4/RTX treated/','EDF11,EDF12/rep1/')]
        hashes=[next(v['sha256'] for k,v in members[n].items() if k in ('unrar','unar','bsdtar') and 'sha256' in v) for n in names]
        assert hashes[0]==hashes[1];alias.append(dict(round=r,members=names,sha256=hashes[0]))
    moment_paths={'npc':'stationarity-diagnostic-sep15/real_split_moments.npz',
        'edf10rep1':'edf10-replication-sep15/moments_edf10.npz',
        'edf10rep2':'replication-cd20-sep15/moments_edf10rep2.npz',
        'fig4untreated':'replication-cd20-sep15/moments_fig4untreated.npz',
        **{f'rtx_rep{r}':f'replication-edf1112-sep15/moments_rep{r}.npz' for r in (1,2,3)},
        'sdc_raw':'independent-field-sep15/moments_raw.npz','sdc_corrected':'independent-field-sep15/moments_corrected.npz'}
    for row in rows:
        tag=row['tag'];v=read(str((P/f'characterization-{tag}.json').relative_to(ROOT)))
        reference_arrays=np.load(P/f'reference-{tag}.npz')
        current=reference_arrays['full_code'];old=np.load(ROOT/'outputs'/moment_paths[tag])['pooled']
        assert current.shape==old.shape
        numeric=dict(max_abs_difference=float(abs(current-old).max()),
                    relative_l2=float(np.linalg.norm(current-old)/np.linalg.norm(old)),
                    passed=bool(np.allclose(current,old,rtol=1e-10,atol=1e-7)))
        if tag in [r['tag'] for r in cd20_registry()] or tag.startswith('sdc_'):
            legacy=replay_legacy_sdc(row) if tag.startswith('sdc_') else replay_legacy_code(row)
            replay=dict(max_abs_difference=float(abs(legacy-old).max()),
                        relative_l2=float(np.linalg.norm(legacy-old)/np.linalg.norm(old)),
                        passed=bool(np.allclose(legacy,old,rtol=1e-10,atol=1e-7)))
            np.savez_compressed(P/f'legacy-replay-{tag}.npz',pooled=legacy)
            v['float64_centered_phase_sensitivity']=numeric
        else:replay=numeric
        v['full_code_moment_replay']=replay;v['spatial_intervals']=spatial_intervals(tag)
        if tag=='npc':
            v['historical_code_variants']=read('outputs/stationarity-diagnostic-sep15/real_stationarity.json')['variants']
        elif tag.startswith('sdc_'):
            v['historical_code_variants']=row['old']['regimes'][tag.replace('sdc_','')]
        else:v['historical_code_variants']=row['old']['code_split']
        empirical=np.load(P/f'empirical-{tag}.npz');events=empirical['events']
        candidates=[];scatter=[];nearest=[]
        for lo,hi in zip(empirical['offsets'][:-1],empirical['offsets'][1:]):
            e=events[lo:hi];labels=e[:,6].astype(int);centroid=[];sizes=[]
            for label in np.unique(labels):
                m=labels==label;pan=e[m,3].astype(int);c=np.bincount(pan,minlength=4)
                candidates.append(c);res=e[m,:2].copy();df=0
                for t in np.unique(pan):
                    k=pan==t;res[k]-=res[k].mean(0);df+=max(int(k.sum())-1,0)
                if df>=4:scatter.append(float((res**2).sum()/(2*df)))
                centroid.append(e[m,:2].mean(0));sizes.append(int(m.sum()))
            centers=np.asarray(centroid)[np.asarray(sizes)>=5]
            if len(centers)>1:
                from scipy.spatial import cKDTree
                nearest.extend(cKDTree(centers).query(centers,k=2)[0][:,1])
        v['empirical_conditions']=dict(candidate_clusters=len(candidates),
            per_panel_event_count_quantiles=np.quantile(candidates,[.1,.5,.9],axis=0),
            eligible_per_panel_event_count_quantiles=np.quantile(reference_arrays['support'][:,2:6],[.1,.5,.9],axis=0),
            eligible_per_panel_localization_count_quantiles=np.quantile(reference_arrays['support'][:,6:10],[.1,.5,.9],axis=0),
            reported_event_variance_quantiles=np.quantile(events[:,4],[.1,.5,.9]),
            event_localization_multiplicity_quantiles=np.quantile(events[:,5],[.1,.5,.9]),
            within_panel_per_axis_scatter_quantiles=np.quantile(scatter,[.1,.5,.9]),
            candidate_neighbor_distance_quantiles=np.quantile(nearest,[.1,.5,.9]),
            scatter_note='Inferred clusters before source QC; scatter can include source merges and within-panel motion. It is not a measured independent-event noise law.')
        if tag in [r['tag'] for r in cd20_registry()]:v['response_sensitivity']=responses(tag)
        out.append(v)
    save(P/'comparison.json',dict(distinct_acquisitions=8,views=9,alias_verification=alias,acquisitions=out))
    with (P/'comparison.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['view','sources','localization_coverage','event_coverage','oracle_unweighted_raw','oracle_event_weight_raw',
            'oracle_localization_weight_raw','full_code_legacy_T','selected_code_legacy_T','diagonal_legacy_T','historical_source_replay','historical_moment_replay'])
        for a in out:w.writerow([a['tag'],a['passing_sources'],a['localization_coverage'],a['event_coverage'],
            *[a['estimates'][k]['raw'] for k in ('unweighted','event_count','localization_count')],
            a['full_code']['total_arithmetic'],a['selected_code']['total_arithmetic'],a['assigned_diagonal']['total_arithmetic'],
            a['reference_replay']['passed'] if a['reference_replay'] else 'no historic raw oracle',a['full_code_moment_replay']['passed']])
    print((P/'comparison.csv').read_text())

if __name__=='__main__':main()
