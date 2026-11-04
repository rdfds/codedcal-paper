"""Readable reports from retained workflow records; no adapter execution."""
import csv
import json
from collections import Counter
from pathlib import Path

from .workflow_decisions import bank_decision
from .workflow_validation import write_json


def write_report(output):
    output=Path(output)
    gate=json.loads((output/'qualification.json').read_text())
    plan=json.loads((output/'plan.json').read_text())
    trial_root=output/'trials';records=[];policies={};native_rows=[]
    active=plan.get('selected_strategy','mixed_conventional')
    def scoped(decision):
        decision=dict(decision)
        decision['decision']={'accept':'all_specified_passed','reject':'observed_failure','unsupported':'unsupported'}[decision['decision']]
        return decision
    if trial_root.exists():
        for path in sorted(trial_root.glob('*/native.json')):
            native=json.loads(path.read_text());policy=path.parent.name;trials=[]
            for p in sorted(path.parent.glob('*.json')):
                if p.name=='native.json':continue
                trial=json.loads(p.read_text());trials.append(trial);records.append(trial)
            by_id={r['challenge']['id']:r for r in trials}
            def decision(ids):
                return scoped(bank_decision([by_id.get(i,dict(assessment=dict(status='unrun',failed=None))) for i in ids]))
            policies[policy]=dict(native_workflow_qualified=native['native_workflow_qualified'],
                native_nonblank_calls=native.get('nonblank_calls'),native_coverage=native.get('coverage'),
                completed_pairs=len(trials),heldout=decision(plan['heldout']),fractional=decision(plan['fractional']),
                strategies={k:decision(v) for k,v in plan['strategies'].items()},
                specified_challenges=decision(plan['strategies'].get(active,[])))
            if not native['native_workflow_qualified']:
                policies[policy]['heldout']=dict(decision='unsupported',reason='native_workflow_not_qualified')
                policies[policy]['specified_challenges']=dict(decision='unsupported',reason='native_workflow_not_qualified')
            native_rows.append(native)
    state_path=trial_root/'state.json'
    state=json.loads(state_path.read_text()) if state_path.exists() else None
    rows=[]
    for call in gate['calls']:
        d=call.get('diagnostics',{})
        classes=call.get('catalog',{}).get('classes')
        blank_count=None if classes is None else sum(k in gate.get('blank_classes',[]) for k in classes)
        rows.append(dict(policy=gate.get('native_policy','identity'),role='control:'+call['role'],
            status=call['status'],reason=call.get('error'),candidate_coverage=call.get('coverage'),
            predicted_calls=None if classes is None else len(classes)-blank_count,blank_calls=blank_count,
            blank_proxy=d.get('blank_proxy'),wall_seconds=call.get('wall_seconds')))
    for native in native_rows:
        rows.append(dict(policy=native['policy'],role='policy_native',status=native['status'],
            reason=native.get('error') or (None if native['native_workflow_qualified'] else 'native_workflow_not_qualified'),
            candidate_coverage=native.get('coverage'),predicted_calls=native.get('nonblank_calls'),
            blank_calls=native.get('diagnostics',{}).get('blank_calls'),
            blank_proxy=native.get('diagnostics',{}).get('blank_proxy'),wall_seconds=native.get('wall_seconds')))
    for record in records:
        a=record['assessment'];candidate,reference=record['runs']
        row=dict(policy=record['policy'],role=record['challenge']['role'],challenge=record['challenge']['id'],
            label=record['challenge']['label'],status=a['status'],reason=a.get('reason'),failed=a.get('failed'),
            candidate_coverage=candidate.get('coverage'),reference_coverage=reference.get('coverage'),
            candidate_error=candidate.get('error'),reference_error=reference.get('error'),
            blank_proxy=candidate.get('diagnostics',{}).get('blank_proxy'),
            oracle_native_identity_loss=record.get('oracle_vs_native',{}).get('identity_loss'),
            oracle_native_composition_tv=record.get('oracle_vs_native',{}).get('composition_tv'),
            actual_workflow_executions=record['actual_workflow_executions'],
            wall_seconds=sum(r.get('wall_seconds',0.) for r in record['runs']))
        for key in ['reference_calls','predicted_calls','blank_calls','identity_loss','composition_tv','yield_ratio','spatial_first_identity_loss','spatial_first_supported']:
            row[key]=a.get(key)
        rows.append(row)
    fields=['policy','role','challenge','label','status','reason','failed','reference_calls','predicted_calls',
        'identity_loss','composition_tv','yield_ratio','candidate_coverage','reference_coverage','blank_calls',
        'blank_proxy','spatial_first_identity_loss','spatial_first_supported','oracle_native_identity_loss',
        'oracle_native_composition_tv','candidate_error','reference_error','actual_workflow_executions','wall_seconds']
    with (output/'recovery.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    execution_records=list(gate['calls'])+native_rows+[r for trial in records for r in trial['runs']]
    def accounting(calls):
        return dict(recorded_workflow_executions=sum(r.get('actual_workflow_executions',r.get('diagnostics',{}).get('actual_workflow_executions',1)) for r in calls),
            exact_result_cache_hits=sum(bool(r.get('diagnostics',{}).get('exact_result_cache_hit')) for r in calls),
            measured_union_seconds=sum(r.get('wall_seconds',0.) for r in calls),
            calibration_executions_reported=sum(r.get('diagnostics',{}).get('actual_calibration_executions',0) for r in calls),
            decoder_executions_reported=sum(r.get('diagnostics',{}).get('actual_decoder_executions',0) for r in calls),
            calls_without_execution_diagnostics=sum('actual_workflow_executions' not in r.get('diagnostics',{}) for r in calls))
    summary=dict(stage='controls_failed' if not gate['passed'] else (state['stage'] if state else 'native_controls_only'),
        selected_strategy=active,design_audit=plan.get('design_audit'),
        outcome_scope='Only the specified digital challenges; all_specified_passed is not a general adequacy or safety conclusion',
        native_policy=gate.get('native_policy','identity'),native_controls_passed=gate['passed'],
        fractional_input_supported=gate['fractional_input_supported'],policies=policies,
        logical_workflow_dispatches=gate['logical_workflow_executions']+(state['logical_workflow_dispatches'] if state else 0),
        accounting=accounting(execution_records),
        per_policy_measured_union_runtime={p:accounting([n for n in native_rows if n['policy']==p]+[r for t in records if t['policy']==p for r in t['runs']]) for p in policies},
        assessment_statuses=dict(Counter(r['assessment']['status'] for r in records)),
        scope='Finite controlled digital challenges, not physical accuracy, RNA truth or unknown-error risk. Union runtimes are not separately measured strategy costs.')
    write_json(output/'summary.json',summary)
    lines=['Workflow challenge report', 'Status: '+summary['stage'],
        'Native policy: '+summary['native_policy'], 'Native controls passed: '+str(gate['passed']),
        'Native nonblank calls: '+str(gate['native_nonblank_calls']),
        'Fractional input supported: '+str(gate['fractional_input_supported'])]
    for p,data in policies.items():
        lines.append(p+': native gate '+str(data['native_workflow_qualified'])+'; specified challenges '+data['specified_challenges']['decision']+'; completed pairs '+str(data['completed_pairs']))
        lines.extend('  '+k+': '+v['decision'] for k,v in data['strategies'].items())
        if plan['heldout']:lines.append('  Additional specified bank: '+data['heldout']['decision'])
    lines.extend(['Logical dispatched queries: '+str(summary['logical_workflow_dispatches']),
        'Recorded workflow executions: '+str(summary['accounting']['recorded_workflow_executions']),summary['outcome_scope'],summary['scope'],
        'See recovery.csv for recovery, composition, support, execution failures and unavailable outcomes.'])
    (output/'summary.txt').write_text('\n'.join(lines)+'\n')
    return summary
