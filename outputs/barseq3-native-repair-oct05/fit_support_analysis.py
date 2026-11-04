"""Interpret both fitting domains with all failures and all original cases."""
import numpy as np
import recovery_analysis as A
from condition import HERE,read,write,sha

rows=[];fields={};complete=True
primary=read(HERE/'recovery/summary.json')
for i in [0,1]:
    folder=HERE/'fit-support'/f'field-{i}';fields[str(i)]={}
    for decoder in A.DECODERS:
        cases=[]
        for kind in A.KINDS:
            gp=folder/f'{kind}-geometry.json';geo=read(gp) if gp.exists() else {}
            values={};quality={}
            for policy in A.POLICIES:
                p=folder/f'{kind}__{policy}__{decoder}.json';r=read(p) if p.exists() else dict(status='missing')
                rows.append(dict(field=i,kind=kind,policy=policy,decoder=decoder,status=r['status'],quality=r.get('quality')))
                if r['status']!='completed':complete=False;continue
                assert r['fit_support_sha256']==sha(HERE/'fit_support.py')
                assert r['specification_sha256']==sha(HERE/'FIT_SUPPORT_SENSITIVITY.md')
                values[policy]={rad:m['gene_class_identity_f1'] for rad,m in r['metrics'].items()};quality[policy]=r['quality']
            cases.append(dict(kind=kind,feature=geo.get('feature'),recommend=geo.get('recommend'),scores=values,quality=quality))
        record=dict(cases=cases)
        if all(len(c['scores'])==5 and c['recommend'] in c['scores'] for c in cases):
            means={p:float(np.mean([c['scores'][p]['5'] for c in cases[1:]])) for p in A.POLICIES}
            means['selected']=float(np.mean([c['scores'][c['recommend']]['5'] for c in cases[1:]]))
            record.update(mean_nonzero_identity_f1=means,
                harmful_recommendations=sum(c['scores'][c['recommend']]['5']-c['scores']['common_iterative']['5']<=-.05 for c in cases),
                decision_diagnostics=A.confusion([c['scores']['population_template']['5']-c['scores']['common_iterative']['5']>=.05 for c in cases],[c['recommend']=='population_template' for c in cases]),
                primary_selected_f1=primary['fields'][str(i)]['decoders'][decoder]['mean_nonzero_identity_f1']['selected'],
                proxy_target_misses={p:sum(not c['quality'][p]['qualified'] for c in cases) for p in A.POLICIES})
        else:complete=False
        fields[str(i)][decoder]=record
assert len(rows)==60
counts={s:sum(r['status']==s for r in rows) for s in ['completed','failed','missing']}
complete=complete and counts['completed']==60
write(HERE/'fit-support/summary.json',dict(counts=counts,planned=60,full_cohort_supported=complete,fields=fields,
    scope='Post-hoc complete-cohort interior-fit sensitivity; same preparations and algorithm reference, not new biological confirmation',analysis_sha256=sha(__file__)))
write(HERE/'fit-support/all-arms.json',rows)
lines=['# Interior-fit sensitivity','',f"Completed {counts['completed']}/60; failed {counts['failed']}; missing {counts['missing']}.",
    f'Full-cohort comparison supported: {complete}.','',
    '| Field | Decoder | Original rule F1 | Interior rule F1 | Interior strong round F1 | Interior author model F1 | Harmful selections |',
    '|---|---|---:|---:|---:|---:|---:|']
for i,decoders in fields.items():
    for d,r in decoders.items():
        if 'mean_nonzero_identity_f1' not in r:continue
        m=r['mean_nonzero_identity_f1'];lines.append(f"| {i} | {d} | {r['primary_selected_f1']:.4f} | {m['selected']:.4f} | {m['common_iterative']:.4f} | {m['barcode_model']:.4f} | {r['harmful_recommendations']}/3 |")
lines+=['','Means use both nonzero cases; zero cases remain in decision/harm diagnostics and full ledgers. Fitting uses only observed interior pixels, with template support checked on every iteration. Smaller fit domains change information as well as padding exposure. Original results remain unchanged. Byte-identical decoder inputs can reuse primary calls, with provenance recorded; they are not independent solver replications. Neither cohort supplies independently known RNA identities or measures temporal T. All three radii, sequence identities, composition and point-proxy/support misses remain per arm.']
(HERE/'fit-support/RESULTS.md').write_text('\n'.join(lines)+'\n')
print(counts,'full cohort',complete,flush=True)
