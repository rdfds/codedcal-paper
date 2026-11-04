"""Replay retained catalogs under the original reference-first eligibility rule.

Never refits or rewrites raw rows. It fixes the documented v1 ordering defect.
"""
from pathlib import Path
import json,copy
import numpy as np
from codedcal.workflow_validation import Catalog,WorkflowResult,assess,compare,canonical


def adjudicate_trial(row,blank_classes,minimum_reference=10):
    result=copy.deepcopy(row);a,b=result['runs'];old=copy.deepcopy(row['assessment'])
    reference=None
    if b['status']=='completed':reference=Catalog.from_dict(b['catalog'])
    nref=0 if reference is None else int(np.sum(~np.isin(reference.classes,blank_classes)))
    if b['status']!='completed':evaluation=dict(status='unsupported',reason='oracle_reference_execution_failed',failed=None)
    elif b['coverage']<.99:evaluation=dict(status='unsupported',reason='oracle_reference_support_failed',failed=None)
    elif nref<minimum_reference:evaluation=dict(status='unsupported',reason='insufficient_reference_calls',failed=None,reference_calls=nref)
    elif a['status']!='completed':evaluation=dict(status='workflow_failed',reason='challenged_execution_failed',failed=True,reference_calls=nref)
    else:
        prediction=WorkflowResult(Catalog.from_dict(a['catalog']),np.array(a['corrections']),a['coverage'],a['diagnostics'])
        evaluation=assess(reference,prediction,blank_classes,minimum_reference=minimum_reference)
        evaluation['matching_sensitivity']={str(rad):compare(reference,prediction.catalog,blank_classes,rad) for rad in [1.5,5.]}
        residual=canonical(np.array(a['corrections'])+np.array(result['challenge']['maps'])-np.array(b['corrections']))
        evaluation['known_probe_map_rms']=float(np.sqrt(np.mean(np.sum(residual**2,axis=1))))
    for metric in [evaluation]+list(evaluation.get('matching_sensitivity',{}).values()):
        if metric.get('matching',{}).get('unresolved_components',0):
            metric['spatial_first_identity_loss']=None;metric['spatial_first_supported']=False
        elif 'spatial_first_identity_loss' in metric:metric['spatial_first_supported']=True
    result['assessment']=evaluation;result['raw_v1_assessment']=old
    result['eligibility_adjudication_changed']=(old.get('failed'),old.get('status'))!=(evaluation.get('failed'),evaluation.get('status'))
    return result
