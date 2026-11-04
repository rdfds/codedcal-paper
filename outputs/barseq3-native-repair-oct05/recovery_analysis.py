"""All-cohort interpretation; incomplete reference gates never create transfer claims."""
from pathlib import Path
import numpy as np
from condition import HERE,read,write,sha

POLICIES=['none','common_iterative','population_template','barcode_model','oracle']
DECODERS=['correlation','regression'];KINDS=['zero','common','color']

def confusion(truth,pred):
    return dict(tp=sum(a and b for a,b in zip(truth,pred)),tn=sum(not a and not b for a,b in zip(truth,pred)),
        fp=sum(not a and b for a,b in zip(truth,pred)),fn=sum(a and not b for a,b in zip(truth,pred)))

def main():
    results={};all_arms=[];complete=True
    for index in [0,1]:
        folder=HERE/'recovery'/f'field-{index}';completion=folder/'completion.json'
        if not completion.exists():
            results[str(index)]=dict(stage='not_completed');complete=False;continue
        state=read(completion)
        if not state.get('reference_qualified',False):
            results[str(index)]=state;complete=False;continue
        field={}
        for decoder in DECODERS:
            cases=[]
            for kind in KINDS:
                geo=read(folder/(kind+'-geometry.json'));scores={};quality={}
                for policy in POLICIES:
                    p=folder/f'{kind}__{policy}__{decoder}.json'
                    if not p.exists():complete=False;all_arms.append(dict(index=index,kind=kind,policy=policy,decoder=decoder,status='missing'));continue
                    arm=read(p);assert arm['recovery_sha256']==sha(HERE/'recovery.py')
                    all_arms.append(dict(index=index,kind=kind,policy=policy,decoder=decoder,status=arm['status'],quality=arm.get('quality')))
                    if arm['status']!='completed':complete=False;continue
                    scores[policy]={r:m['gene_class_identity_f1'] for r,m in arm['metrics'].items()};quality[policy]=arm['quality']
                cases.append(dict(kind=kind,feature=geo['feature'],recommend=geo['recommend'],scores=scores,quality=quality,
                    complete=len(scores)==len(POLICIES)))
            record=dict(cases=cases)
            if all(c['complete'] for c in cases):
                record['mean_nonzero_identity_f1']={p:float(np.mean([c['scores'][p]['5'] for c in cases if c['kind']!='zero'])) for p in POLICIES}
                record['mean_nonzero_identity_f1']['selected']=float(np.mean([c['scores'][c['recommend']]['5'] for c in cases if c['kind']!='zero']))
                # Include zero in decision diagnostics, retaining unnecessary/harmful recalibration.
                record['decision_diagnostics']={r:confusion([c['scores']['population_template'][r]-c['scores']['common_iterative'][r]>=.05 for c in cases],[c['recommend']=='population_template' for c in cases]) for r in ['1.5','3','5']}
                record['harmful_recommendations']=sum(c['scores'][c['recommend']]['5']-c['scores']['common_iterative']['5']<=-.05 for c in cases)
                record['zero_selected_identity_f1']=cases[0]['scores'][cases[0]['recommend']]['5']
                record['error_grouping_effect']={p:cases[1]['scores'][p]['5']-cases[2]['scores'][p]['5'] for p in POLICIES}
                record['proxy_target_misses']={p:sum(not c['quality'][p]['qualified'] for c in cases) for p in POLICIES}
            field[decoder]=record
        results[str(index)]=dict(stage=state['stage'],decoders=field)
    # Add unattempted entries explicitly to the 60-arm ledger.
    keys={(a['index'],a['kind'],a['policy'],a['decoder']) for a in all_arms}
    for i in [0,1]:
        for k in KINDS:
            for p in POLICIES:
                for d in DECODERS:
                    if (i,k,p,d) not in keys:all_arms.append(dict(index=i,kind=k,policy=p,decoder=d,status='not_run'))
    assert len(all_arms)==60
    counts={s:sum(a['status']==s for a in all_arms) for s in ['completed','failed','missing','not_run']}
    complete=complete and counts['completed']==60
    summary=dict(planned_arms=60,counts=counts,full_cohort_comparison_supported=complete,fields=results,
        biological_scope='Two fields on one BARseq3 preparation; separate from older CSHL benchmark. Adapter development/spatial holdout, not prospective intervention.',
        endpoint='Consistency of native algorithm gene-class outputs after controlled integer displacement; not independent RNA truth',
        perturbation='Exactly identical full offset multiset and covariance, grouped by round or color; rank-one covariance, seed505',
        analysis_sha256=sha(__file__))
    write(HERE/'recovery/summary.json',summary);write(HERE/'recovery/all-arms.json',all_arms)
    lines=['# Public matched-distribution recovery','',summary['endpoint']+'.','',
        f"Completed {counts['completed']}/60 planned arms; failed {counts['failed']}, missing {counts['missing']}, not run {counts['not_run']}.",
        '**Full-cohort comparison supported: '+str(complete)+'.**','',
        '| Field | Decoder | Strong common F1 | Frozen recommendation F1 | Author model F1 | Harmful recommendations (including zero) |',
        '|---|---|---:|---:|---:|---:|']
    for i,field in results.items():
        for d,r in field.get('decoders',{}).items():
            if 'mean_nonzero_identity_f1' in r:
                m=r['mean_nonzero_identity_f1'];lines.append(f"| {i} | {d} | {m['common_iterative']:.4f} | {m['selected']:.4f} | {m['barcode_model']:.4f} | {r['harmful_recommendations']} / 3 |")
    lines+=['','Both nonzero patterns have the exact same complete marginal displacement distribution: twenty +(4,2) and twenty -(4,2). Only grouping by round versus color changes. All failures and blank/support target misses remain in all-arms.json; decision diagnostics include the zero case and three matching radii. More conditions on these fields would not add biological preparations. No true FDR, randomized temporal validation or guaranteed journal acceptance follows. The original annotated study supplies a different, incomplete author-annotation endpoint and must be interpreted separately.']
    (HERE/'recovery/RESULTS.md').write_text('\n'.join(lines)+'\n');print(counts,'full cohort:',complete,flush=True)

if __name__=='__main__':main()
