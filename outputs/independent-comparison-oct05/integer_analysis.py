"""Interpret every planned integer-control arm; no refitting or case selection."""
import json
from pathlib import Path
import numpy as np
import analyze as study

HERE=Path(__file__).resolve().parent
POLICIES=study.POLICIES

def main():
    counts=dict(planned=120,completed=0,failed_or_missing=0);results={};rows=[]
    for split in ['development','validation']:
        folder=HERE/'integer-control'/split
        paths=[folder/(f'{kind}-a{dose}.json') for kind,dose in [('zero',0),('round',2),('round',4),('color',2),('color',4)]]
        assert all(p.exists() for p in paths),'Incomplete condition cohort'
        cases=[study.read(p) for p in paths]
        zero=cases[0]
        assert all(c['integer_control_sha256']==study.sha(HERE/'integer_control.py') for c in cases)
        results[split]={}
        for decoder in study.DECODERS:
            base=zero['arms']['none__'+decoder]
            for case in cases:
                oracle=case['arms'].get('oracle__'+decoder)
                assert oracle is not None,'Oracle failed'
                assert oracle['xy']==base['xy'] and oracle['words']==base['words'],'Oracle changes calls despite exact integer interior'
            valid=[]
            for case in cases:
                values={}
                for policy in POLICIES:
                    key=policy+'__'+decoder
                    arm=case['arms'].get(key)
                    counts['completed']+=arm is not None;counts['failed_or_missing']+=arm is None
                    rows.append(dict(split=split,condition=case['condition'],decoder=decoder,policy=policy,
                        status='completed' if arm is not None else 'failed_or_missing',
                        identity_f1=arm['metrics']['5']['identity_f1'] if arm else None))
                    if arm is not None:values[policy]=arm['metrics']['5']['identity_f1']
                if not case['condition'].startswith('zero'):
                    valid.append(dict(condition=case['condition'],feature=case['feature'],recommend=case['recommend'],
                        identity_f1=values,complete=all(p in values for p in POLICIES)))
            complete=[c for c in valid if c['complete']]
            assert len(complete)==len(valid),'Full-comparison interpretation withheld: missing arm'
            gains=[c['identity_f1']['population_template']-c['identity_f1']['common_iterative'] for c in complete]
            predicted=[c['recommend']=='population_template' for c in complete]
            means={p:float(np.mean([c['identity_f1'][p] for c in complete])) for p in POLICIES}
            means['selected']=float(np.mean([c['identity_f1'][c['recommend']] for c in complete]))
            results[split][decoder]=dict(cases=valid,oracle_calls_equal_zero=True,
                mean_identity_f1=means,confusion=study.confusion([g>=.05 for g in gains],predicted),
                harmful_recommendations=sum(c['identity_f1'][c['recommend']]-c['identity_f1']['common_iterative']<=-.05 for c in complete))
    out=dict(arms=counts,results=results,scope='Post-hoc exposed-crop sensitivity with integer warps and a common prepared decoder window; not independent confirmation')
    study.write(HERE/'integer-control/summary.json',out)
    study.write(HERE/'integer-control/all-arms.json',rows)
    lines=['# Integer-only resampling sensitivity','',out['scope']+'.','',
        f"Completed {counts['completed']}/120 planned arms; failed/missing {counts['failed_or_missing']}.",
        'Oracle calls exactly reproduce zero/no-correction calls in both families and both crops.','',
        '| Crop | Decoder | Strong round-only F1 | Frozen recommendation F1 | Author model F1 | Harmful recommendations |',
        '|---|---|---:|---:|---:|---:|']
    for split in results:
        for decoder,r in results[split].items():
            m=r['mean_identity_f1'];lines.append(f"| {split} | {decoder} | {m['common_iterative']:.4f} | {m['selected']:.4f} | {m['barcode_model']:.4f} | {r['harmful_recommendations']} / 4 |")
    lines+=['','Integer application removes interpolation blur; identical oracle calls also control decoder-domain/gain/threshold changes. Applied fitted maps are quantized, so this is a stress test rather than a claim about subpixel performance. The original 312 arms, false recommendation, resampling caveat and external qualification failure remain. All three matching-radius metrics and complete maps/calls are retained per condition.']
    (HERE/'integer-control/RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({s:{d:r['mean_identity_f1'] for d,r in v.items()} for s,v in results.items()},indent=2))

if __name__=='__main__':main()
