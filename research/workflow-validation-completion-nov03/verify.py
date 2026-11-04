"""Verify frozen completion and retained evidence without fitting any workflow."""
from pathlib import Path
import argparse
import fcntl
import json
import sys

H = Path(__file__).resolve().parent
R = H.parents[1]
P = R/'research/workflow-validation-controller-nov02'
sys.path.insert(0, str(H))
from prepare import verify_manifest
from recover import sha, missing_pairs


def main(seal=False):
    f = json.loads((H/'freeze.json').read_text())
    original = json.loads((P/'freeze.json').read_text())
    assert sha(P/'freeze.json') == f['original_freeze_sha256']
    allowed_changes = {'inherited_call_count', 'total_reserved', 'continuation_of', 'wall_seconds_limit',
        'minimum_start_free_bytes', 'source_comparative_outcomes_examined', 'comparative_outcomes_examined', 'preceding_commit'}
    for key, value in original.items():
        if key not in allowed_changes:
            assert f[key] == value, key
    assert f['original_reservation_preserved'] == original['total_reserved'] == 1150
    assert f['total_reserved'] == 1151 == 1150+f['additional_interrupted_retry_reserved'] <= f['maximum_total'] == 1600
    assert f['maximum_new_dispatches'] == 6 and len(f['missing_pairs']) == 3
    for pins in [f['code_sha256'], f['recovery_code_sha256']]:
        for name, digest in pins.items():
            assert sha(R/name) == digest, name
    assert sha(R/f['input']['path']) == f['input']['sha256']
    prefix = json.loads((H/'prefix-manifest.json').read_text())
    for row in prefix['files']:
        assert sha(R/row['path']) == row['sha256'] == sha(R/row['source'])
    assert verify_manifest(P/'closure-manifest.json') == 216
    prior_counts = {}
    for name in prefix['previous_sealed_files_verified']:
        path = R/'research'/name
        prior_counts[name] = verify_manifest(path)
    source_analysis = (P/'analyze.py').read_text()
    assert sha(P/'analyze.py') == f['analysis_source_sha256']
    for a, b in f['analysis_text_only_replacements'].items():
        source_analysis = source_analysis.replace(a, b)
    assert (H/'analyze.py').read_text() == source_analysis

    report = dict(passed=True, prior_manifests_verified=prior_counts,
        inherited_checkpoint_files_verified=len(prefix['files']), original_scientific_design_unchanged=True,
        recovery_reservation=1151, no_new_data_or_thresholds=True, physical_validation_added=False)
    state_path = H/'run/state.json'
    if state_path.exists():
        lock = (H/'run/collector.lock').open('r')
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        s = json.loads(state_path.read_text())
        assert s['stage'] in {'completed', 'stopped_at_limit', 'stopped_with_error'}
        assert s['freeze_sha256'] == sha(H/'freeze.json')
        assert s['logical_call_starts'] <= 64
        calls = [json.loads(p.read_text()) for p in sorted((H/'run/calls').glob('*.json'))]
        assert [c['call'] for c in calls] == list(range(s['logical_call_starts']))
        old_failed = calls[57]
        assert old_failed['stage'] == 'failed' and sha(H/'run/calls/0057.json') == f['retry']['sha256']
        if len(calls) > 58:
            assert calls[58]['input_fingerprint'] == old_failed['input_fingerprint']
            assert calls[58]['retry_of'] == f['retry']
        for index, c in enumerate(calls[58:]):
            assert c['challenge'] == f['missing_pairs'][index//2]['challenge']['id']
            assert c['policy'] == 'group_search'
            if index % 2 and c['stage'] == 'completed':
                assert c['diagnostics']['exact_result_cache_hit'] and c['diagnostics']['actual_workflow_executions'] == 0
        complete_calls = [c for c in calls if c['stage'] == 'completed']
        assert len(complete_calls) == s['complete_calls']
        actual = sum(c['diagnostics']['actual_workflow_executions'] for c in complete_calls)
        hits = sum(c['diagnostics']['exact_result_cache_hit'] for c in complete_calls)
        assert not list((H/'run').rglob('input.npz'))
        remaining = missing_pairs(f, H/'run')
        summary = json.loads((H/'summary.json').read_text())
        assert len(summary['missing']) == len(remaining)
        assert summary['complete_integer_matrix'] == (s['stage']=='completed' and not remaining)
        assert summary['adequate_cases'] == 0 and summary['failed_cases'] == 2
        assert not summary['general_utility_supported'] and not summary['source_panel_strict_increment']
        assert not summary['comparison_rates_estimable']
        if s['stage'] == 'completed':
            assert not remaining and s['completed_pairs'] == 24 and s['completed_new_pairs'] == 3
            assert len(calls) == 64
        report.update(worker_stage=s['stage'], completed_pairs=s['completed_pairs'], missing_pairs=len(remaining),
            full_chain_dispatches=len(calls), completed_actual_workflows=actual, exact_cache_hits=hits,
            failed_or_partial_dispatches=len(calls)-len(complete_calls), new_dispatches=len(calls)-58,
            complete_integer_matrix=summary['complete_integer_matrix'], adequate_cases=0, failed_cases=2,
            comparative_benefit_established=False, general_utility_supported=False)
        manifest = H/'completion-manifest.json'
        if seal:
            assert not manifest.exists(), 'Never overwrite a sealed completion'
            paths = [p for p in H.rglob('*') if p.is_file() and '__pycache__' not in p.parts
                     and p.name not in {'completion-manifest.json', 'verification.json'}]
            rows = [dict(path=str(p.relative_to(R)), bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(paths)]
            manifest.write_text(json.dumps(dict(files=rows), indent=2)+'\n')
        if manifest.exists():
            report['new_sealed_files_verified'] = verify_manifest(manifest)
    else:
        assert not seal
        report.update(stage='frozen_not_started', new_actual_workflows=0)
    (H/'verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--seal', action='store_true')
    main(parser.parse_args().seal)
