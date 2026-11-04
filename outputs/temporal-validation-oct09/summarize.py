"""Recount frozen outputs, cluster paired code draws by physical field."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import xlogy
import validation as V
from common import save, sha


def numeric(rows,key):
    return np.array([np.nan if r.get(key) is None else r[key] for r in rows],float)


def cluster_interval(matrix, rng, valid=None):
    ids=rng.integers(len(matrix),size=(2000,len(matrix)))
    if valid is None: draws=matrix[ids].mean(axis=(1,2))
    else:
        den=valid[ids].sum(axis=(1,2)); num=np.where(valid,matrix,0)[ids].sum(axis=(1,2))
        draws=np.divide(num,den,out=np.full(len(ids),np.nan),where=den>0)
    finite=draws[np.isfinite(draws)]
    return np.quantile(finite,[.025,.975]).tolist() if len(finite) else [None,None]


def bounded_cluster_interval(mean, fields, alpha=.05):
    """Post-hoc sensitivity: finite-sample Chernoff bounds for iid bounded field means.

    Within-field code assignments may have arbitrary dependence. This never
    replaces the frozen percentile-bootstrap gate.
    """
    level=np.log(2/alpha)/fields
    def kl(p):
        return xlogy(mean,mean/p)+xlogy(1-mean,(1-mean)/(1-p))
    if mean==0:return [0.,float(1-np.exp(-level))]
    if mean==1:return [float(np.exp(-level)),1.]
    lower=brentq(lambda p:kl(p)-level,np.finfo(float).tiny,mean)
    upper=brentq(lambda p:kl(p)-level,mean,1-np.finfo(float).eps)
    return [float(lower),float(upper)]


def statistics(rows, rep_count, assignments, label):
    assert len(rows)==rep_count*assignments
    rows=sorted(rows,key=lambda r:(r['rep'],r['assignment']))
    ratio=numeric(rows,'ratio'); usable=np.array([r['usable'] for r in rows]); covers=np.array([r['covers'] for r in rows])
    available=np.array([r['positive_reference'] for r in rows]) & np.isfinite(ratio)
    physical=rows[0]['physical_target']; observable=rows[0]['observable_target']
    valid=available.reshape(rep_count,assignments)
    rng=np.random.default_rng(V.seed(label,'bootstrap'))
    joint=(usable&covers).reshape(rep_count,assignments)
    ci=cluster_interval(joint,rng)
    conditional=cluster_interval(covers.reshape(rep_count,assignments),rng,usable.reshape(rep_count,assignments))
    bias=float(np.mean(ratio[available])-physical) if available.any() else None
    observable_bias=float(np.mean(ratio[available])-observable) if available.any() else None
    width=numeric(rows,'upper')-numeric(rows,'lower')
    d_valid=np.array([r['D']['valid'] for r in rows]);d_quiet=np.array([r['D']['quiet'] for r in rows])
    bad=available&(abs(ratio-physical)>.01)
    nominal_failure='none' if ci[0] is not None and ci[0]<=.95<=ci[1] else ('above_nominal' if ci[0] is not None and ci[0]>.95 else 'below_nominal')
    rejects_zero=usable&(numeric(rows,'upper')<1)
    qualified=bool(bias is not None and abs(bias)<=.01 and usable.mean()>=.90 and joint.mean()>=.90
                   and ci[0] is not None and ci[0]<=.95<=ci[1])
    stat=dict(records=len(rows),physical_fields=rep_count,assignments_per_field=assignments,
              physical_target=physical,observable_target=observable,
              mean_ratio=float(np.mean(ratio[available])) if available.any() else None,
              ratio_bias=bias,observable_ratio_bias=observable_bias,
              ratio_bias_interval=cluster_interval((ratio-physical).reshape(rep_count,assignments),rng,valid),
              ratio_rmse=float(np.sqrt(np.mean((ratio[available]-physical)**2))) if available.any() else None,
              positive_reference_availability=float(available.mean()),
              available_ratio_quantiles=np.quantile(ratio[available],[0,.01,.05,.5,.95,.99,1]) if available.any() else [],
              usable_fraction=float(usable.mean()),set_types=dict(Counter(r['kind'] for r in rows)),
              unconditional_coverage=float(covers.mean()),
              unconditional_coverage_interval=cluster_interval(covers.reshape(rep_count,assignments),rng),
              usable_joint_coverage=float(joint.mean()),usable_joint_coverage_interval=ci,
              coverage_bootstrap_nominal_failure_direction=nominal_failure,
              coverage_bounded_cluster_interval=bounded_cluster_interval(float(joint.mean()),rep_count),
              finite_interval_excludes_zero_motion_fraction=float(rejects_zero.mean()),
              usable_covers_and_excludes_zero_fraction=float((rejects_zero&covers).mean()),
              conditional_usable_coverage=float(covers[usable].mean()) if usable.any() else None,
              conditional_usable_coverage_interval=conditional,
              median_usable_interval_width=float(np.median(width[usable])) if usable.any() else None,
              D_invalid=int((~d_valid).sum()),D_alerts=int((d_valid&~d_quiet).sum()),
              D_alert_fraction_all=float((d_valid&~d_quiet).mean()),
              D_alert_fraction_valid=float((~d_quiet[d_valid]).mean()) if d_valid.any() else None,
              bad_physical_errors=int(bad.sum()),bad_errors_quiet_D=int((bad&d_quiet).sum()),
              bad_error_quiet_fraction=float((bad&d_quiet).sum()/bad.sum()) if bad.any() else None,
              qualified=qualified,
              failed_gates=[k for k,ok in dict(ratio_bias=bias is not None and abs(bias)<=.01,
                                              usable_fraction=usable.mean()>=.90,
                                              usable_coverage=joint.mean()>=.90,
                                              nominal_coverage_in_bootstrap=ci[0] is not None and ci[0]<=.95<=ci[1]).items() if not ok])
    return stat


def summarize(require_complete=True):
    frozen=json.loads((V.P/'freeze.json').read_text())
    assert frozen['code_sha256']==sha(V.P/'validation.py')
    records=[];worlds=[];artifacts={};field_total=0;group_total=0
    for tag in V.TAGS:
        for model in V.MODELS:
            for world in V.WORLDS:
                folder=V.P/'worlds'/f'{tag}-{model}-{world}'
                statepath=folder/'state.json'
                if not statepath.exists():
                    if require_complete: raise RuntimeError(f'Missing {folder}')
                    continue
                state=json.loads(statepath.read_text())
                if require_complete: assert state['complete'] and state['completed']==state['expected'],folder
                assert state['freeze_sha256']==sha(V.P/'freeze.json')
                artifacts[str(statepath.relative_to(V.P))]=sha(statepath)
                truthfile=folder/'truth.npz'; artifacts[str(truthfile.relative_to(V.P))]=sha(truthfile)
                target_array=np.load(truthfile)['population']
                data=[];clusters=[];seconds=[]
                for entry in state['replicates']:
                    rep=entry['rep'];path=folder/entry['filename'];recordpath=folder/f'field-{rep:03d}.json'
                    assert sha(path)==entry['sha256'];assert sha(recordpath)==entry['record_sha256']
                    artifacts[str(path.relative_to(V.P))]=sha(path);artifacts[str(recordpath.relative_to(V.P))]=sha(recordpath)
                    r=json.loads(recordpath.read_text());data.extend(r['cases']);clusters.extend(r['clustering']);seconds.append(r['seconds'])
                    group_total+=2*max(json.loads((V.P/f'catalog-{tag}.json').read_text())['groups'],512)
                field_total+=state['completed']
                summaries=[]
                for scope in ('observed','standardized'):
                    for method,ncode in (('coded',2),('known_source',1)):
                        for q in V.Q:
                            rows=[r for r in data if r['scope']==scope and r['method']==method and r['q']==float(q)]
                            stat=statistics(rows,state['completed'],ncode,f'{tag}-{model}-{world}-{scope}-{method}-{q}')
                            stat.update(tag=tag,model=model,world=world,scope=scope,method=method,q=float(q),
                                        primary=world in V.PRIMARY,complete=state['complete'],group_count=rows[0]['group_count'])
                            iq=int(np.flatnonzero(V.Q==q)[0])
                            for ti,tn in ((0,'baseline'),(1,'observable')):
                                perwindow=target_array[ti,iq,:,0]/target_array[ti,iq,:,1]
                                stat[f'known_source_D_population_{tn}']=float(perwindow[0]/perwindow[1])
                            stat['common_motion_D_null_declared']=world in ('zero','homogeneous','group_common')
                            summaries.append(stat);records.append(stat)
                cluster_summary=None
                if clusters:
                    countkeys=('candidate_clusters','passing_clusters','impure_passing','events','localizations',
                               'passing_localizations','pure_kernel_mass','total_kernel_mass')
                    cluster_summary={key:sum(c[key] for c in clusters) for key in countkeys}
                    cluster_summary['draws']=len(clusters)
                    cluster_summary['localization_coverage']=cluster_summary['passing_localizations']/cluster_summary['localizations']
                    cluster_summary['pure_kernel_mass_coverage']=cluster_summary['pure_kernel_mass']/cluster_summary['total_kernel_mass']
                    cluster_summary['impure_passing_fraction']=cluster_summary['impure_passing']/cluster_summary['passing_clusters'] if cluster_summary['passing_clusters'] else None
                    cluster_summary['rejections']=dict(sum((Counter(c['rejections']) for c in clusters),Counter()))
                worlds.append(dict(tag=tag,model=model,world=world,complete=state['complete'],physical_fields=state['completed'],
                                   runtime_seconds=sum(seconds),clustering=cluster_summary))
                print('summarized',tag,model,world,flush=True)
    for tag in V.TAGS:
        for extension in ('json','npz'):
            path=V.P/f'catalog-{tag}.{extension}'
            if path.exists():artifacts[path.name]=sha(path)
    for name in ('PROTOCOL.md','IMPLEMENTATION.md','UNCERTAINTY.md','controls.json','controls.py','validation.py','freeze.json','summarize.py','execution-state.json','execution.log','verification-attempt-01.log'):
        artifacts[name]=sha(V.P/name)
    execution=json.loads((V.P/'execution-state.json').read_text())
    summary=dict(execution=execution,settings_frozen_date='2026-10-09',complete=len(worlds)==126 and all(w['complete'] for w in worlds),
                 distinct_acquisitions=8,views=9,worlds=worlds,records=records,
                 physical_fields=field_total,coded_fields=field_total*2,coded_group_draws=group_total,
                 conditional_only=True,physical_validation=False,physical_fields_semantics='independent simulated noise+motion realizations, not acquired images',
                 code_sha256=sha(__file__),freeze_sha256=sha(V.P/'freeze.json'))
    save(V.P/'summary.json',summary)
    simplekeys=['tag','model','world','primary','scope','method','q','group_count','physical_fields','records',
                'physical_target','observable_target','mean_ratio','ratio_bias','ratio_rmse',
                'positive_reference_availability','usable_fraction','unconditional_coverage','usable_joint_coverage',
                'conditional_usable_coverage','median_usable_interval_width','finite_interval_excludes_zero_motion_fraction','usable_covers_and_excludes_zero_fraction','coverage_bootstrap_nominal_failure_direction','D_invalid','D_alerts',
                'bad_physical_errors','bad_errors_quiet_D','known_source_D_population_baseline','known_source_D_population_observable','common_motion_D_null_declared','qualified']
    with (V.P/'summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=simplekeys);w.writeheader();w.writerows({k:r.get(k) for k in simplekeys} for r in records)
    save(V.P/'artifact-manifest.json',dict(artifacts=artifacts,complete=summary['complete'],physical_fields=field_total,coded_group_draws=group_total))
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');a=p.parse_args();summarize(not a.partial)
