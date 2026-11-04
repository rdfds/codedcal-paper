"""Analyze a stopped frozen run without refitting or changing thresholds."""
from pathlib import Path
import collections
import csv
import fcntl
import json
import sys
import numpy as np
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path[:0]=[str(H),str(R/'tool/src')]
from codedcal.workflow_decisions import bank_decision,confusion
from codedcal.workflow_validation import write_json


def main():
    f=json.loads((H/'freeze.json').read_text());out=H/'run'
    lock=(out/'collector.lock').open('r');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state=json.loads((out/'state.json').read_text())
    supervisor=H/'supervisor-state.json'
    assert state['stage'] in {'completed','stopped_at_limit','stopped_with_error'} or (supervisor.exists() and json.loads(supervisor.read_text())['stage']=='stopped'),'Do not analyze a running collection'
    rows=[];missing=[];profiles=[];names=list(f['strategies'])
    for policy in f['policies']:
        q=out/'qualification'/(policy+'.json');gate=json.loads(q.read_text()) if q.exists() else None
        row=dict(policy=policy,preparation=f['preparation'],native_eligible=bool(gate and gate['primary_native_eligible']),native_quality=gate['quality'] if gate else None)
        if gate is None:missing.append(policy+'/qualification')
        if not row['native_eligible']:
            row.update(stage='native_excluded' if gate else 'unrun',truth=dict(decision='unsupported',reason='native_gate_excluded' if gate else 'unrun_native'))
            for name in names:row[name]=dict(decision='unsupported',reason=row['truth']['reason'])
            rows.append(row);continue
        trials={}
        for ch in f['challenges']:
            path=out/'trials'/policy/(ch['id']+'.json')
            if not path.exists():missing.append(policy+'/'+ch['id']);continue
            trial=json.loads(path.read_text())
            if ch['id'] in f['fractional']:
                assert trial['status']=='not_run' and trial['reason']==f['fractional_exclusion']
                profiles.append(dict(policy=policy,role='fractional',challenge=ch['id'],label=ch['label'],status='not_run',reason=trial['reason']))
                continue
            assert trial['challenge']==ch and trial['policy']==policy
            trials[ch['id']]=trial;a=trial['assessment']
            profiles.append(dict(policy=policy,role=ch['role'],challenge=ch['id'],label=ch['label'],status=a['status'],reason=a.get('reason'),
                failed=a.get('failed'),identity_loss=a.get('identity_loss'),composition_tv=a.get('composition_tv'),
                reference_calls=a.get('reference_calls'),predicted_calls=a.get('predicted_calls'),yield_ratio=a.get('yield_ratio'),
                candidate_coverage=trial['runs'][0].get('coverage'),reference_coverage=trial['runs'][1].get('coverage'),
                blank_proxy=trial['runs'][0].get('diagnostics',{}).get('blank_proxy'),integer_oracle_native_exact=trial.get('integer_oracle_native_exact')))
        def decision(ids):return bank_decision([trials.get(i,dict(assessment=dict(failed=None,status='unrun'))) for i in ids])
        row.update(stage='measured',truth=decision(f['heldout']),**{name:decision(ids) for name,ids in f['strategies'].items()});rows.append(row)
    complete=state['stage']=='completed' and not missing
    conf={name:confusion(rows,name) for name in names};known=[r for r in rows if r['truth']['decision']!='unsupported']
    adequate=sum(r['truth']['decision']=='accept' for r in known);failed=len(known)-adequate
    a=conf['mechanism_weighted'];b=conf['mixed_conventional']
    comparison_estimable=adequate>0 and failed>0
    increment=bool(complete and comparison_estimable and a['adequate_accepted']>0 and
        a['failed_detection_rate']>b['failed_detection_rate'] and a['false_alert_rate']<=b['false_alert_rate'] and a['decision_coverage']>=b['decision_coverage'])
    runtime=collections.defaultdict(lambda:dict(measured_union_seconds=0.,actual_workflows=0,actual_calibrations=0,actual_decoders=0,exact_cache_hits=0,failed_or_partial_calls=0))
    call_count=0
    for path in sorted((out/'calls').glob('*.json')):
        call=json.loads(path.read_text());call_count+=1;v=runtime[call['policy']]
        v['measured_union_seconds']+=call.get('wall_seconds',0.)
        if call['stage']!='completed':v['failed_or_partial_calls']+=1;continue
        d=call['diagnostics'];v['actual_workflows']+=d['actual_workflow_executions'];v['actual_calibrations']+=d['actual_calibration_executions'];v['actual_decoders']+=d['actual_decoder_executions'];v['exact_cache_hits']+=int(d['exact_result_cache_hit'])
    result=dict(complete=complete,complete_integer_matrix=complete,worker_stage=state['stage'],rows=rows,missing=missing,confusion=conf,
        supported_heldout_cases=len(known),adequate_cases=adequate,failed_cases=failed,comparison_rates_estimable=comparison_estimable,
        source_panel_strict_increment=increment,
        general_utility_supported=bool(increment and adequate>=f['utility_minimum_adequate'] and failed>=f['utility_minimum_failed'] and len({r['preparation'] for r in known})>=f['utility_minimum_preparations']),
        independent_preparations=1,per_policy_measured_query_union=dict(runtime),new_logical_calls_started=state['logical_call_starts']-f['inherited_call_count'],inherited_logical_dispatches=f['inherited_call_count'],full_chain_logical_dispatches=state['logical_call_starts'],
        recorded_call_checkpoints=call_count,conservative_total_reservation=f['total_reserved'],
        inherited_native_controls_are_new_executions=False,fractional_panel_eligible=False,
        scope='One retrospectively selected/native-qualified public preparation, frozen comparative digital challenges. Native catalog preservation is not independent molecular truth or physical validation.')
    write_json(H/'summary.json',result)
    with (H/'decisions.csv').open('w',newline='') as stream:
        fields=['policy','native_eligible','native_nonblank_calls','native_blank_calls','native_blank_proxy_upper95','heldout','mixed_conventional','mechanism_weighted']
        w=csv.DictWriter(stream,fieldnames=fields);w.writeheader()
        for r in rows:
            q=r['native_quality'] or {};w.writerow(dict(policy=r['policy'],native_eligible=r['native_eligible'],native_nonblank_calls=q.get('nonblank_calls'),native_blank_calls=q.get('blank_calls'),native_blank_proxy_upper95=q.get('blank_proxy_upper95'),heldout=r['truth']['decision'],**{k:r[k]['decision'] for k in names}))
    fields=list(dict.fromkeys(k for p in profiles for k in p))
    if fields:
        with (H/'profiles.csv').open('w',newline='') as stream:
            w=csv.DictWriter(stream,fieldnames=fields);w.writeheader();w.writerows(profiles)
    lines=['# Completion of the frozen public-source qualification check','',
        'Completed integer matrix: '+str(complete)+'. Supported held-out cases: '+str(len(known))+'. Adequate: '+str(adequate)+'. Failed: '+str(failed)+'.',
        '', '| Suite | Failed detected | Adequate accepted | False alerts | Unsupported decisions |', '| --- | ---: | ---: | ---: | ---: |']
    for name,c in conf.items():lines.append('| '+name+' | '+str(c['failed_rejected'])+' | '+str(c['adequate_accepted'])+' | '+str(c['adequate_rejected'])+' | '+str(c['unsupported_decision'])+' |')
    lines += ['', 'Strict descriptive increment on this source panel: '+str(increment)+'. General utility supported: '+str(result['general_utility_supported'])+'.',
        '', 'This is one additional preparation with four assigned calibration policies. Methods, challenges and pixels are not independent acquisitions. Missing adequate or failed examples prevents a useful rate comparison. A tie, reject-all or incomplete matrix is not superiority.',
        '', 'Every assigned policy and native-quality exclusion is in decisions.csv. Fractional challenge outcomes were not run: the retained composed native control fails the unchanged .10 proxy-quality requirement. Native references remain algorithmic references, not molecular truth.',
        '', 'Budget reservation: '+str(f['total_reserved'])+'/1600. Full-chain logical calls including all 58 inherited dispatches and the preserved interruption: '+str(state['logical_call_starts'])+'. Actual executions/cache hits and per-policy measured query unions are in summary.json. All inherited checkpoints and imported native controls are not newly measured runs; strategy runtimes are not separately measured.',
        '', 'No new source requests, paid jobs, reserved data, physical-validation claim, release or submission. All old negative findings remain unchanged.']
    (H/'RESULTS.md').write_text('\n'.join(lines)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['rows','per_policy_measured_query_union','missing']},indent=2))


if __name__=='__main__':main()
