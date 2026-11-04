"""Adjudicate stable retained decoder outputs without another prediction."""
import csv
import fcntl
from replay_common import *
from codedcal.workflow_decisions import bank_decision,confusion


def main():
    f=load(H/'freeze.json');out=H/'run'
    lock=(out/'collector.lock').open('r');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state=load(out/'state.json');assert state['stage'] in {'completed','stopped_at_limit','stopped_with_error'}
    rows=[];profiles=[];missing=[];runtime={p:{m:0. for m in ['historical','anchored']} for p in f['policies']}
    for e in f['cases']:
        for mode in ['historical','anchored']:
            path=out/mode/(e['id']+'.json')
            if not path.exists():missing.append(mode+'/'+e['id']);continue
            r=load(path);runtime[e['policy']][mode]+=r.get('seconds',0.)
            if r['stage']!='completed':missing.append(mode+'/'+e['id']+':'+r['stage'])
    complete=state['stage']=='completed' and not missing and state['complete_predictions']==52
    for policy in f['policies']:
        oldgate=load(P/'run/qualification'/(policy+'.json'))
        modes={}
        for mode in ['historical','anchored']:
            trials={}
            for e in f['cases']:
                if e['policy']!=policy:continue
                path=out/mode/(e['id']+'.json')
                r=load(path) if path.exists() else dict(stage='unrun')
                if e['challenge'] is None:
                    native=r;continue
                trials[e['challenge']['id']]=r
                a=r.get('assessment',{})
                profiles.append(dict(policy=policy,mode=mode,challenge=e['challenge']['id'],label=e['challenge']['label'],
                    role=e['challenge']['role'],stage=r['stage'],status=a.get('status'),failed=a.get('failed'),
                    identity_loss=a.get('identity_loss'),composition_tv=a.get('composition_tv'),
                    predicted_calls=a.get('predicted_calls'),reference_calls=a.get('reference_calls'),coverage=r.get('coverage'),
                    blank_proxy_upper95=r.get('quality',{}).get('blank_proxy_upper95')))
            def decide(ids):
                return bank_decision([dict(assessment=trials.get(i,{}).get('assessment',dict(status='unrun',failed=None))) for i in ids])
            q=native.get('quality',{});a=native.get('assessment',{})
            eligible=bool(native['stage']=='completed' and q.get('nonblank_calls',0)>=50 and
                q.get('blank_proxy_upper95') is not None and q['blank_proxy_upper95']<=.1 and native.get('coverage',0)>=.99)
            compatible=a.get('failed') is False
            modes[mode]=dict(native_quality=q,native_eligible=eligible,native_fixed_reference_compatible=compatible,
                native_assessment=a,heldout=decide(f['heldout']),
                **{name:decide(ids) for name,ids in f['strategies'].items()},
                integer_passes=sum(r.get('assessment',{}).get('failed') is False for r in trials.values()),
                integer_failures=sum(r.get('assessment',{}).get('failed') is True for r in trials.values()),
                integer_unsupported=sum(r.get('assessment',{}).get('failed') is None for r in trials.values()))
        rescued=bool(complete and oldgate['primary_native_eligible'] and modes['anchored']['native_eligible'] and
            modes['anchored']['native_fixed_reference_compatible'] and modes['historical']['heldout']['decision']=='reject' and
            modes['anchored']['heldout']['decision']=='accept')
        rows.append(dict(policy=policy,preparation=f['preparation'],original_native_eligible=oldgate['primary_native_eligible'],
            historical=modes['historical'],anchored=modes['anchored'],development_workflow_rescued=rescued))
    changes=[]
    for p in profiles:
        if p['mode']!='anchored':continue
        before=next(r for r in profiles if r['policy']==p['policy'] and r['challenge']==p['challenge'] and r['mode']=='historical')
        changes.append(dict(policy=p['policy'],challenge=p['challenge'],role=p['role'],before_failed=before['failed'],after_failed=p['failed'],
            identity_loss_change=None if before['identity_loss'] is None or p['identity_loss'] is None else p['identity_loss']-before['identity_loss']))
    conf={}
    for mode in ['historical','anchored']:
        decisions=[dict(truth=r[mode]['heldout'],**{name:r[mode][name] for name in f['strategies']}) for r in rows]
        conf[mode]={name:confusion(decisions,name) for name in f['strategies']}
    result=dict(complete=complete,worker_stage=state['stage'],missing_or_failed=missing,policies=rows,case_changes=changes,
        development_workflows_rescued=sum(r['development_workflow_rescued'] for r in rows),
        failures_to_passes=sum(c['before_failed'] is True and c['after_failed'] is False for c in changes),
        passes_to_failures=sum(c['before_failed'] is False and c['after_failed'] is True for c in changes),
        strategy_confusion=conf,measured_per_policy_union_seconds=runtime,decoder_attempts=state['decoder_attempts'],
        actual_decoder_predictions_started=state['decoder_predictions_started'],completed_predictions=state['complete_predictions'],
        fixed_codebook_index_builds=state['fixed_codebook_index_builds'],exact_cache_hits=0,
        registration_optimizations=0,gain_optimizations=0,software_decoder_predictions=4,
        total_conservative_reservation=f['total_conservative_reservation'],maximum_total=f['maximum_total'],
        elapsed_wall_seconds=state['elapsed_wall_seconds'],elapsed_monotonic_seconds=state['elapsed_monotonic_seconds'],
        peak_rss_bytes=state['peak_rss_bytes'],independent_preparations=1,independent_confirmation=False,
        general_utility_supported=False,physical_validation_added=False,
        fractional_panel_excluded=True,original_native_exclusions=['round','frame_template'],
        minimum_general_utility_unchanged=dict(adequate=5,failed=5,preparations=2),
        scope='Post-hoc causal intervention only on final common gain scale, with saved geometry/relative gains and original native references fixed.')
    target=H/'summary.json'
    if target.exists():assert load(target)==result,'Offline adjudication changed'
    else:write_json(target,result)
    with (H/'profiles.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(profiles[0]) if profiles else ['stage']);writer.writeheader();writer.writerows(profiles)
    lines=['# Fixed-fit gain-scale replay','',f'Complete replay: {complete}. Development workflows rescued: {result["development_workflows_rescued"]}/2.',
        '', '| Workflow | Historical integer passes | Anchored integer passes | Historical held-out | Anchored held-out | Native compatible |',
        '| --- | ---: | ---: | --- | --- | --- |']
    for r in rows:
        a,b=r['historical'],r['anchored']
        lines.append(f'| {r["policy"]} | {a["integer_passes"]}/12 | {b["integer_passes"]}/12 | {a["heldout"]["decision"]} | {b["heldout"]["decision"]} | {b["native_fixed_reference_compatible"]} |')
    lines+=['',f'Individual failures rescued: {result["failures_to_passes"]}; passing cases regressed: {result["passes_to_failures"]}.',
        '',f'Executions: {state["decoder_predictions_started"]}/52 decoder predictions started; {state["complete_predictions"]} completed, plus four tiny software controls. No registration/gain optimization or result-cache substitution. One fixed codebook index was constructed. Conservative reservation1215/1600.',
        '', 'All historical reproductions must pass before any anchored case. Catalogs, all24 assessments, blank proxies, support, matching-radius sensitivities and preserved exclusions are retained. Original native references and .05 loss gates are unchanged.',
        '', 'Both conventional and mechanism-weighted suites remain visible in summary.json; no per-strategy runtime was separately measured. Decoder-only execution counts do not imply cost equality with full workflows.',
        '', 'This is exposed development on one preparation. It cannot establish independent confirmation, useful qualification across preparations, RNA truth or physical accuracy. Five adequate/five failed workflows across two preparations are still required for the original utility claim. Fractional challenges and the two originally excluded native workflows remain excluded.']
    if complete and not result['development_workflows_rescued']:
        lines+=['','The proposed common-scale intervention did not rescue either workflow. Close this gain-scale-benefit branch. Keep the explicit gain-units contract without claiming a recovery improvement. No additional experiment is started.']
    (H/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:result[k] for k in ['complete','development_workflows_rescued','failures_to_passes','passes_to_failures','actual_decoder_predictions_started','general_utility_supported']},indent=2))


if __name__=='__main__':main()
