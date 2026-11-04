"""Read retained outcomes without rerunning/altering methods or thresholds."""
from pathlib import Path
import sys,json,hashlib
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path.insert(0,str(R/'tool/src'))
from codedcal.workflow_validation import write_json
from codedcal.workflow_decisions import bank_decision,confusion,leave_preparation_out,qc_union_lopo
freeze=json.loads((H/'development-freeze.json').read_text());OUT=H/'development';rows=[];missing=[];runs=[]
for e in freeze['cases']:
 for decoder in e['decoders']:
  for policy in freeze['policies']:
   path=OUT/e['name']/decoder/policy;nativepath=path/'native.json'
   if not nativepath.exists():missing.append([e['name'],decoder,policy]);continue
   native=json.loads(nativepath.read_text());trials={}
   for ch in e['plan']['challenges']:
    p=path/(ch['id']+'.json')
    if p.exists():
     trial=json.loads(p.read_text());trials[ch['id']]=trial;runs+=trial['runs']
   def get(ids):return [trials.get(k,dict(assessment=dict(failed=None,status='unrun'))) for k in ids]
   r=dict(name=e['name'],preparation=e['preparation'],decoder=decoder,policy=policy,native=native,
    truth=bank_decision(get(e['plan']['heldout'])),fractional=bank_decision(get(e['plan']['fractional'])))
   for key,ids in e['plan']['strategies'].items():
    r[key]=bank_decision(get(ids));r[key+'_geometric']=r[key]['geometric_score'];r[key+'_yield']=r[key]['yield_score']
   diag=native.get('diagnostics',{});r['native_blank_proxy']=diag.get('blank_proxy');r['native_fit_loss']=diag.get('fit_loss');r['native_inverse_yield']=1/max(diag.get('nonblank_calls',0),1) if diag else None
   # Native failure/insufficient usable catalog prevents trivial full-workflow qualification.
   if native['status']!='completed' or native.get('coverage',0)<.99 or diag.get('nonblank_calls',0)<50:
    for key in e['plan']['strategies']:r[key]['decision']='reject';r[key]['native_gate_failed']=True
    r['truth']['decision']='unsupported';r['truth']['native_gate_failed']=True
   rows.append(r)
for feature in ['native_blank_proxy','native_fit_loss','native_inverse_yield']+[k+'_'+s for k in ['conventional','mixed_conventional','balanced_random','readout_aware'] for s in ['geometric','yield']]:
 predictions=leave_preparation_out(rows,feature)
 for r,p in zip(rows,predictions):r[feature+'_lopo']=p
for r,p in zip(rows,qc_union_lopo(rows)):r['combined_native_qc_lopo']=p
keys=['combined_native_qc_lopo','conventional','mixed_conventional','balanced_random','readout_aware']+['native_blank_proxy_lopo','native_fit_loss_lopo','native_inverse_yield_lopo']+[k+'_'+s+'_lopo' for k in ['conventional','mixed_conventional','balanced_random','readout_aware'] for s in ['geometric','yield']]
summary=dict(rows=rows,missing_workflows=missing,complete=len(missing)==0 and json.loads((OUT/'state.json').read_text())['stage']=='completed',
 confusion={k:confusion(rows,k) for k in keys},per_preparation={p:{k:confusion([r for r in rows if r['preparation']==p],k) for k in keys} for p in sorted({r['preparation'] for r in rows})},
 measured_actual_trial_workflow_executions=sum(r.get('actual_workflow_executions',1) for r in runs),
 measured_trial_wall_seconds=sum(r.get('wall_seconds',0) for r in runs),calibration_executions=sum(r.get('diagnostics',{}).get('actual_calibration_executions',0) for r in runs),decoder_executions=sum(r.get('diagnostics',{}).get('actual_decoder_executions',0) for r in runs),
 scope='Finite digital challenges, descriptive prep-level comparisons; unsupported rules are not successful qualifications; cache hits carry lookup time, not inherited runtime')
write_json(H/'development-summary.json',summary)
print('COMPLETE',summary['complete'],'rows',len(rows),'missing',len(missing))
for k in keys:print(k,summary['confusion'][k])
