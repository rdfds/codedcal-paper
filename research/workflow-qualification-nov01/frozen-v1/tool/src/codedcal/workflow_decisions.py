"""Finite challenge-bank decisions; no natural-error or biological truth claims."""
import numpy as np
import itertools


def bank_decision(rows,tolerance=.05):
    """A valid failure suffices to reject; acceptance needs all supported passes."""
    assessments=[r['assessment'] for r in rows]
    failures=[r for r in assessments if r.get('failed') is True]
    supported=[r for r in assessments if r.get('failed') is not None]
    values=[max(r.get('identity_loss',1. if r.get('failed') else 0.),r.get('composition_tv',0.)) for r in supported]
    state='reject' if failures else ('accept' if len(supported)==len(rows) and rows else 'unsupported')
    return dict(decision=state,score=max(values) if values else None,supported_challenges=len(supported),planned_challenges=len(rows),
        failures=len(failures),tolerance=tolerance,geometric_score=max((r.get('known_probe_map_rms',0.) for r in supported),default=None),
        yield_score=max((abs(r.get('yield_ratio',1.)-1.) for r in supported),default=None))


def confusion(records,prediction_key,truth_key='truth'):
    counts={k:0 for k in ['failed_rejected','failed_accepted','adequate_rejected','adequate_accepted','unsupported_decision','unsupported_truth']}
    for r in records:
        truth=r[truth_key]['decision'];pred=r[prediction_key]['decision']
        if truth=='unsupported':counts['unsupported_truth']+=1;continue
        if pred=='unsupported':counts['unsupported_decision']+=1;continue
        key=('failed_' if truth=='reject' else 'adequate_')+('rejected' if pred=='reject' else 'accepted');counts[key]+=1
    adequate=counts['adequate_accepted']+counts['adequate_rejected'];failed=counts['failed_rejected']+counts['failed_accepted']
    supported=sum(counts[k] for k in ['failed_rejected','failed_accepted','adequate_rejected','adequate_accepted'])
    counts.update(false_alert_rate=counts['adequate_rejected']/adequate if adequate else None,
        failed_detection_rate=counts['failed_rejected']/failed if failed else None,adequate_acceptance_rate=counts['adequate_accepted']/adequate if adequate else None,
        decision_coverage=supported/(supported+counts['unsupported_decision']) if supported+counts['unsupported_decision'] else None)
    return counts


def fit_threshold(training,score_key,max_false_alert=.1,minimum_adequate=5):
    """Train one observable scalar threshold; count missing scores as abstentions.

    Empirical training constraints are not a population confidence guarantee.
    This intentionally abstains if exposed adequate examples are insufficient.
    """
    usable=[r for r in training if r['truth']['decision']!='unsupported' and r.get(score_key) is not None and np.isfinite(r[score_key])]
    adequate=[r for r in usable if r['truth']['decision']=='accept']
    if len(adequate)<minimum_adequate:return dict(supported=False,reason='fewer_than_five_adequate_training_workflows',adequate_training=len(adequate))
    candidates=[np.nextafter(max(r[score_key] for r in usable),np.inf)]+sorted(set(r[score_key] for r in usable))
    valid=[]
    for cut in candidates:
        fp=sum(r[score_key]>cut for r in adequate)/len(adequate)
        bad=[r for r in usable if r['truth']['decision']=='reject'];tp=sum(r[score_key]>cut for r in bad)/len(bad) if bad else 0.
        if fp<=max_false_alert:valid.append((-tp,fp,float(cut)))
    _,fp,cut=min(valid)
    return dict(supported=True,threshold=cut,training_false_alert_rate=fp,adequate_training=len(adequate),training_workflows=len(usable),population_error_guarantee=False)


def leave_preparation_out(records,score_key):
    predictions=[]
    for row in records:
        training=[r for r in records if r['preparation']!=row['preparation']]
        rule=fit_threshold(training,score_key);score=row.get(score_key)
        decision='unsupported' if not rule['supported'] or score is None or not np.isfinite(score) else ('reject' if score>rule['threshold'] else 'accept')
        predictions.append(dict(decision=decision,rule=rule,score=score,test_preparation_excluded_from_training=True))
    return predictions


def fit_qc_union(training,features=('native_blank_proxy','native_fit_loss','native_inverse_yield')):
    """Strong, small-data QC comparator: train an OR of up to two scalar gates.

    The selected feature family must have >=80% coverage of training workflows.
    Missing blank measurements never become zero-valued evidence of safety.
    """
    known=[r for r in training if r['truth']['decision']!='unsupported'];best=None;best_rule=None
    for n in [1,2]:
        for family in itertools.combinations(features,n):
            usable=[r for r in known if all(r.get(f) is not None and np.isfinite(r[f]) for f in family)]
            adequate=[r for r in usable if r['truth']['decision']=='accept'];bad=[r for r in usable if r['truth']['decision']=='reject']
            if len(adequate)<5 or len(usable)<.8*len(known):continue
            candidates=[sorted({r[f] for r in usable}|{float(np.nextafter(max(r[f] for r in usable),np.inf))}) for f in family]
            for cuts in itertools.product(*candidates):
                def alert(r):return any(r[f]>cut for f,cut in zip(family,cuts))
                fp=sum(map(alert,adequate))/len(adequate)
                if fp>.1:continue
                tp=sum(map(alert,bad))/len(bad) if bad else 0.
                key=(-tp,fp,-len(usable),n,family,cuts)
                if best is None or key<best:
                    best=key;best_rule=dict(supported=True,features=family,thresholds=cuts,training_false_alert_rate=fp,training_failure_detection=tp,adequate_training=len(adequate),training_coverage=len(usable)/len(known),population_error_guarantee=False)
    return best_rule or dict(supported=False,reason='insufficient_adequate_training_or_feature_coverage')


def qc_union_lopo(records):
    result=[]
    for r in records:
        rule=fit_qc_union([t for t in records if t['preparation']!=r['preparation']])
        decision='unsupported'
        if rule['supported'] and all(r.get(f) is not None and np.isfinite(r[f]) for f in rule['features']):
            decision='reject' if any(r[f]>cut for f,cut in zip(rule['features'],rule['thresholds'])) else 'accept'
        result.append(dict(decision=decision,rule=rule,test_preparation_excluded_from_training=True))
    return result
