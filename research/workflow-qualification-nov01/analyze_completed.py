"""Reference-first analysis of a stopped finite run; raw evidence is immutable."""
from pathlib import Path
import sys,json,hashlib,collections,csv
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path[:0]=[str(H),str(R/'tool/src')]
from codedcal.workflow_validation import ImageCase,write_json
from codedcal.workflow_decisions import bank_decision,confusion,leave_preparation_out,qc_union_lopo,fit_threshold,fit_qc_union
from adjudicate import adjudicate_trial
F=json.loads((H/'development-freeze.json').read_text());O=H/'development';state=json.loads((O/'state.json').read_text());assert state['stage']!='running','Never interpret a changing collection as completed evidence'
strategies=['conventional','mixed_conventional','balanced_random','readout_aware'];rows=[];missing=[];changes=[];runtime=collections.defaultdict(lambda:dict(actual_workflows=0,actual_calibrations=0,actual_decoders=0,wall_seconds=0.,calibration_seconds=0.,decoder_seconds=0.,lookup_hits=0));raw_count=0;primary_oracle_mismatches=[];profiles=[]
def charge(policy,r):
 d=r.get('diagnostics',{});v=runtime[policy];v['actual_workflows']+=r.get('actual_workflow_executions',d.get('actual_workflow_executions',1));v['actual_calibrations']+=d.get('actual_calibration_executions',0);v['actual_decoders']+=d.get('actual_decoder_executions',0);v['wall_seconds']+=r.get('wall_seconds',d.get('wall_seconds',0.));v['calibration_seconds']+=d.get('calibration_seconds',0.);v['decoder_seconds']+=d.get('decoder_seconds',0.);v['lookup_hits']+=int(d.get('exact_result_cache_hit',False))
for e in F['cases']:
 p=Path(e['path'])
 if not p.exists():p=H/'development-inputs'/Path(e['path']).name
 assert hashlib.sha256(p.read_bytes()).hexdigest()==e['sha256'];case=ImageCase.load(p)
 for decoder in e['decoders']:
  for policy in F['policies']:
   key=[e['name'],decoder,policy];dest=O/e['name']/decoder/policy;nativepath=dest/'native.json';source_qualified=e['qualified'][decoder]
   row=dict(name=e['name'],preparation=e['preparation'],decoder=decoder,policy=policy,source_qualified=source_qualified,declared_target=case.metadata['target'],support_scope='candidate domain' if 'object_adapter' in case.metadata else 'scored voxels',logical_probe_budget_per_strategy=8)
   if not source_qualified:
    row.update(stage='source_gate_excluded',native_workflow_qualified=False,truth=dict(decision='unsupported',reason='source_gate_excluded'))
    for k in strategies:row[k]=dict(decision='unsupported',reason='source_gate_excluded')
    rows.append(row);continue
   if not nativepath.exists():
    missing.append(key);row.update(stage='unrun',native_workflow_qualified=False,truth=dict(decision='unsupported',reason='native_workflow_unrun'))
    for k in strategies:row[k]=dict(decision='unsupported',reason='native_workflow_unrun')
    rows.append(row);continue
   native=json.loads(nativepath.read_text());row['native']=native;charge(policy,native);trials={}
   for ch in e['plan']['challenges']:
    path=dest/(ch['id']+'.json')
    if not path.exists():missing.append(key+[ch['id']]);continue
    raw=json.loads(path.read_text());assert raw['challenge']==ch and raw['policy']==policy
    trial=adjudicate_trial(raw,case.blank_classes);trials[ch['id']]=trial;raw_count+=1
    a=trial['assessment'];native_floor=trial.get('oracle_vs_native',{})
    profiles.append(dict(input=e['name'],preparation=e['preparation'],decoder=decoder,policy=policy,role=ch['role'],challenge=ch['id'],label=ch['label'],status=a['status'],reason=a.get('reason'),failed=a.get('failed'),reference_calls=a.get('reference_calls'),predicted_calls=a.get('predicted_calls'),identity_loss=a.get('identity_loss'),composition_tv=a.get('composition_tv'),yield_ratio=a.get('yield_ratio'),known_map_rms=a.get('known_probe_map_rms'),candidate_coverage=trial['runs'][0].get('coverage'),reference_coverage=trial['runs'][1].get('coverage'),oracle_native_identity_loss=native_floor.get('identity_loss'),oracle_native_composition_tv=native_floor.get('composition_tv'),actual_workflow_executions=trial['actual_workflow_executions']))
    if trial['eligibility_adjudication_changed']:changes.append(dict(case=key,challenge=ch['id'],old=raw['assessment'],corrected=trial['assessment']))
    for r in raw['runs']:charge(policy,r)
    if ch['role']!='fractional' and native['status']=='completed' and raw['runs'][1]['status']=='completed':
     b=raw['runs'][1]
     if b['catalog']!=native['catalog'] or not np.array_equal(b['corrections'],native['corrections']):primary_oracle_mismatches.append(key+[ch['id']])
   def get(ids):return [trials.get(i,dict(assessment=dict(failed=None,status='unrun'))) for i in ids]
   row.update(stage='measured',truth=bank_decision(get(e['plan']['heldout'])),fractional=bank_decision(get(e['plan']['fractional'])),completed_challenge_pairs=len(trials))
   d=native.get('diagnostics',{});row['native_workflow_qualified']=native['status']=='completed' and native.get('coverage',0)>=.99 and d.get('nonblank_calls',0)>=50
   row['native_blank_proxy']=d.get('blank_proxy');row['native_fit_loss']=d.get('fit_loss');row['native_inverse_yield']=1/max(d.get('nonblank_calls',0),1) if d else None
   for k in strategies:
    row[k]=bank_decision(get(e['plan']['strategies'][k]));row[k+'_geometric']=row[k]['geometric_score'];row[k+'_yield']=row[k]['yield_score']
   if not row['native_workflow_qualified']:
    row['truth']['decision']='unsupported';row['truth']['reason']='native_workflow_not_qualified'
    for k in strategies:row[k]['decision']='reject';row[k]['native_gate_failed']=True
   rows.append(row)
# Preserve every initially planned image-regression exclusion without calling it a failure.
for token in F['known_native_regression_compute_exclusions']:
 name,decoder=token.split('/');prep=next(e['preparation'] for e in F['cases'] if e['name']==name)
 for policy in F['policies']:
  r=dict(name=name,preparation=prep,decoder=decoder,policy=policy,stage='compute_excluded',source_qualified=None,native_workflow_qualified=None,truth=dict(decision='unsupported',reason='compute_excluded'))
  for k in strategies:r[k]=dict(decision='unsupported',reason='compute_excluded')
  rows.append(r)
features=['native_blank_proxy','native_fit_loss','native_inverse_yield']+[k+'_'+m for k in strategies for m in ['geometric','yield']]
for feature in features:
 for r,p in zip(rows,leave_preparation_out(rows,feature)):r[feature+'_lopo']=p
for r,p in zip(rows,qc_union_lopo(rows)):r['combined_native_qc_lopo']=p
methods=strategies+['combined_native_qc_lopo']+[f+'_lopo' for f in features]
known=[r for r in rows if r['truth']['decision']!='unsupported'];adequate=sum(r['truth']['decision']=='accept' for r in known);failed=len(known)-adequate
conf={k:confusion(rows,k) for k in methods};a=conf['readout_aware'];comparators=[conf[k] for k in strategies if k!='readout_aware']
# Only fixed-gate empirical domination, with nontrivial acceptance/coverage, could support an increment.
utility_eligible=state['stage']=='completed' and not missing and adequate>=5 and failed>=5 and len({r['preparation'] for r in known})>=2
strict_improvement=utility_eligible and a['adequate_accepted']>0 and a['failed_detection_rate'] is not None and all(c['failed_detection_rate'] is not None and c['false_alert_rate'] is not None and a['false_alert_rate']<=c['false_alert_rate'] and a['failed_detection_rate']>c['failed_detection_rate'] and a['decision_coverage']>=c['decision_coverage'] for c in comparators)
result=dict(complete=state['stage']=='completed' and not missing,worker_stage=state['stage'],rows=rows,missing=missing,raw_pairs_replayed=raw_count,reference_gate_corrections=changes,integer_oracle_native_mismatches=primary_oracle_mismatches,
 original_planned_workflows=32,compute_feasible_workflows=24,source_qualified_workflows=F['qualified_workflows'],heldout_truth_supported=len(known),adequate_workflows=adequate,failed_workflows=failed,
 empirical_utility_comparison_supported=utility_eligible,strict_increment_over_all_strong_suites=strict_improvement,
 confusion=conf,per_preparation={p:{k:confusion([r for r in rows if r['preparation']==p],k) for k in methods} for p in sorted({r['preparation'] for r in rows})},
 per_policy_measured_union_runtime=dict(runtime),worker_elapsed_monotonic_seconds=state['elapsed_seconds'],logical_execution_reservation=state['logical_charged'],
 logical_workflow_dispatches=2*raw_count+sum('native' in r for r in rows),
 actual_workflow_executions=sum(v['actual_workflows'] for v in runtime.values()),
 actual_calibration_executions=sum(v['actual_calibrations'] for v in runtime.values()),
 actual_decoder_executions=sum(v['actual_decoders'] for v in runtime.values()),
 exact_result_cache_hits=sum(v['lookup_hits'] for v in runtime.values()),
 independent_confirmation_completed=False,physical_validation_added=False,
 scope='Known digital interventions on declared local/conditional catalog targets. Invocation matched, not asserted equal wall-clock seconds. Fields/decoders/challenges are not independent preparations. No natural-risk or RNA-truth claim.')
assert not primary_oracle_mismatches,'Known integer inverse did not reproduce the complete native workflow'
write_json(H/'adjudicated-summary.json',result)
with (H/'workflow-decisions.csv').open('w',newline='') as stream:
 fields=['input','preparation','decoder','policy','stage','source_qualified','native_workflow_qualified','truth']+strategies
 writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
 for r in rows:writer.writerow(dict(input=r['name'],preparation=r['preparation'],decoder=r['decoder'],policy=r['policy'],stage=r['stage'],source_qualified=r['source_qualified'],native_workflow_qualified=r['native_workflow_qualified'],truth=r['truth']['decision'],**{k:r[k]['decision'] for k in strategies}))
if profiles:
 with (H/'recovery-profiles.csv').open('w',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(profiles[0]));writer.writeheader();writer.writerows(profiles)
write_json(H/'frozen-development-rules.json',dict(fixed_output_tolerance=.05,minimum_native=50,minimum_reference=10,minimum_coverage=.99,native_qc_union=fit_qc_union(rows),scalar_rules={f:fit_threshold(rows,f) for f in features},training_preparations=sorted({r['preparation'] for r in known}),confirmation_outcomes_examined=False))
print('Complete:',result['complete'],'supported:',len(known),'adequate:',adequate,'failed:',failed,'strict increment:',strict_improvement)
for k in strategies:print(k,conf[k])
