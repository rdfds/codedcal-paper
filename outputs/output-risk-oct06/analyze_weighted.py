"""Predeclared descriptive risk evaluation; no outcome-selected primary threshold."""
import collections,hashlib,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent/'weighted';ROOT=HERE.parent.parent.parent

def read(p):return json.loads(p.read_text())
def roc(scores,truth):
    data=[(s,y) for s,y in zip(scores,truth) if s is not None];positive=sum(y for s,y in data);negative=len(data)-positive
    if not positive or not negative:return dict(auc=None,sensitivity_at_FPR_01=None,available=len(data),positive=positive,negative=negative)
    auc=sum((a>b)+.5*(a==b) for a,y in data if y for b,z in data if not z)/(positive*negative)
    curve=[]
    for cut in sorted(set(s for s,y in data),reverse=True)+[float('-inf')]:
        tp=sum(y and s>=cut for s,y in data);fp=sum(not y and s>=cut for s,y in data);curve.append(dict(threshold=None if not np.isfinite(cut) else cut,tpr=tp/positive,fpr=fp/negative))
    return dict(auc=auc,sensitivity_at_FPR_01=max([r['tpr'] for r in curve if r['fpr']<=.1]+[0.]),available=len(data),positive=positive,negative=negative,curve=curve,threshold_curve_is_exposed_diagnostic_not_frozen_transfer=True)

def main():
    paths=sorted((HERE/'cases').glob('*/evaluation.json'));assert len(paths)==13,len(paths);records=[];abstained_cases=[]
    for p in paths:
        d=read(p)
        if d['stage']!='complete':abstained_cases.append(dict(case=p.parent.name,reason=d['reason']));continue
        for r in d['rows']:
            risk=r['risk'];actual=r['actual_consistency_error'];records.append(dict(case=d['case'],preparation=d['source'],field=d['field'],kind=d['kind'],decoder=r['decoder'],reference_calls=r['reference_calls'],incumbent_calls=r['incumbent_calls'],truth_failure=r['failed_consistency'],actual=actual,abstained=risk['abstained'],risk=None if risk['abstained'] else max(risk['identity_risk'],risk['composition_risk']),flagged=None if risk['abstained'] else risk['flagged'],ordinary=None if risk['abstained'] else max(risk['ordinary_multistart_identity_disagreement'],risk['ordinary_multistart_composition_disagreement']),full_bank=None if risk['abstained'] else max(risk['unfiltered_bank_identity_disagreement'],risk['unfiltered_bank_composition_disagreement']),fit=r['baseline_fit_loss'],fit_gap=r['baseline_fit_gap'],low_confidence=None if r['baseline_decoder_confidence'] is None else 1-r['baseline_decoder_confidence'],blank_proxy=r['baseline_quality']['blank_proxy'],yield_loss=-r['incumbent_calls'],geometry=None if risk['abstained'] else risk['weighted_geometric_RMS']))
    correlations=[r for r in records if r['decoder']=='correlation'];supported=[r for r in correlations if not r['abstained'] and r['reference_calls']>=10];positive=sum(r['truth_failure'] for r in supported);negative=len(supported)-positive;tp=sum(r['flagged'] and r['truth_failure'] for r in supported);fp=sum(r['flagged'] and not r['truth_failure'] for r in supported);sensitivity=tp/positive if positive else None;fpr=fp/negative if negative else None;coverage=len(supported)/13
    benchmarks={key:roc([r[key] for r in supported],[r['truth_failure'] for r in supported]) for key in ['risk','ordinary','full_bank','fit','fit_gap','low_confidence','blank_proxy','yield_loss','geometry']}
    prep={}
    for p in ['bar','culture','tissue']:
        rows=[r for r in supported if r['preparation']==p];prep[p]=dict(supported=len(rows),failures=sum(r['truth_failure'] for r in rows),misses=sum(r['truth_failure'] and not r['flagged'] for r in rows),false_alerts=sum(not r['truth_failure'] and r['flagged'] for r in rows),field_ids=sorted(set(r['field'] for r in rows)))
    failure_preps=sum(p['failures']>0 for p in prep.values());a=benchmarks['risk']['sensitivity_at_FPR_01'];b=[benchmarks[k]['sensitivity_at_FPR_01'] for k in ['ordinary','full_bank','fit_gap','low_confidence','geometry']];adds_value=a is not None and all(v is not None and a>v for v in b)
    passes=coverage>=.8 and sensitivity is not None and sensitivity>=.9 and fpr is not None and fpr<=.1 and failure_preps>=2 and adds_value
    result=dict(stage='interpreted',posthoc_statistic_repair=True,planned_cases=13,completed_cases=len(paths),planned_first_stage_decoder_cases=20,completed_decoder_cases=len(records),correlation_supported_cases=len(supported),correlation_decision_coverage=coverage,correlation_sensitivity=sensitivity,correlation_false_alert_rate=fpr,failure_preparations=failure_preps,adds_value_over_full_bank_and_ordinary_disagreement=adds_value,first_stage_qualified=passes,nonnegative_extension_authorized_by_scientific_gate=passes,nonnegative_public_propagation_completed=0,independent_preparation_confirmation=False,independent_RNA_truth=False,formal_risk_coverage_guarantee=False,development=True,preparation_counts_are_not_independent_fields=prep,comparators=benchmarks,abstained_cases=abstained_cases,records=records,retained_evaluation_hashes=[dict(file=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]);(HERE/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['records','comparators','retained_evaluation_hashes']},indent=2))
if __name__=='__main__':main()
