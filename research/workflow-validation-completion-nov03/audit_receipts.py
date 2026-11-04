"""Offline native-child, retry and accounting audit; no fits or source access."""
from pathlib import Path
import hashlib
import json

H = Path(__file__).resolve().parent
R = H.parents[1]
P = R/'research/workflow-validation-controller-nov02'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    f = json.loads((H/'freeze.json').read_text())
    summary = json.loads((H/'summary.json').read_text())
    old_call = json.loads((P/'run/calls/0057.json').read_text())
    old_receipt = json.loads((P/'run/children/0022/receipt.json').read_text())
    assert old_call['stage'] == 'failed' and old_receipt['exit_code'] == -15
    assert sha(P/'run/calls/0057.json') == f['retry']['sha256']
    verified = []
    for job, index in enumerate([58, 60, 62]):
        directory = H/'run/children'/f'{job:04}'
        receipt = json.loads((directory/'receipt.json').read_text())
        payload = json.loads((directory/'result.json').read_text())
        call = json.loads((H/'run/calls'/f'{index:04}.json').read_text())
        assert receipt['stage'] == payload['stage'] == call['stage'] == 'completed'
        assert receipt['exit_code'] == 0 and receipt['input_removed']
        assert receipt['local_input_bytes'] <= 64*1024**2 and not (directory/'input.npz').exists()
        assert receipt['result_sha256'] == sha(directory/'result.json')
        assert receipt['input_fingerprint'] == call['input_fingerprint']
        assert call['catalog'] == payload['catalog'] and call['corrections'] == payload['corrections']
        assert call['coverage'] == payload['coverage']
        assert call['diagnostics']['actual_workflow_executions'] == 1
        assert not call['diagnostics']['exact_result_cache_hit']
        assert receipt['child_peak_rss_bytes'] == payload['peak_rss_bytes'] == call['diagnostics']['child_peak_rss_bytes']
        assert json.loads((directory/'config.json').read_text()) == f['adapter_config']
        if job == 0:
            assert receipt['input_sha256'] == old_receipt['input_sha256']
            assert call['input_fingerprint'] == old_call['input_fingerprint']
        reference = json.loads((H/'run/calls'/f'{index+1:04}.json').read_text())
        assert reference['stage'] == 'completed' and reference['diagnostics']['exact_result_cache_hit']
        assert reference['diagnostics']['actual_workflow_executions'] == 0
        assert reference['diagnostics']['child_process_executions'] == 0
        verified.append(dict(job=job, call=index, result_sha256_verified=True,
            input_removed_after_hash_check=True, input_sha256=receipt['input_sha256'],
            input_fingerprint=receipt['input_fingerprint'], child_peak_rss_bytes=receipt['child_peak_rss_bytes']))
    for key, upstream in [('bardensr_sha256', 'bardensr-upstream.nosync'), ('serval_sha256', 'serval-upstream.nosync')]:
        for name, digest in f[key].items():
            assert sha(R.parent/upstream/name) == digest, name
    calls = [json.loads(p.read_text()) for p in sorted((H/'run/calls').glob('*.json'))]
    complete = [c for c in calls if c['stage'] == 'completed']
    assert len(calls) == 64 and len(complete) == 63
    assert len(list((H/'run/children').iterdir())) == 3
    assert sum(c['diagnostics']['actual_workflow_executions'] for c in complete) == 29
    assert sum(c['diagnostics']['exact_result_cache_hit'] for c in complete) == 34
    assert summary['complete_integer_matrix'] and not summary['missing']
    report = dict(passed=True, new_actual_native_workflows=3, new_exact_cache_hits=3, new_dispatches=6,
        all_three_native_results_hash_verified=True, all_native_configurations_identical_to_freeze=True,
        interrupted_input_full_serialization_sha256_identical=True, old_interruption_preserved=True,
        exact_inverse_references_reused_without_fitting=True, upstream_code_pins_unchanged=True, children=verified)
    path = H/'receipt-verification.json'
    if path.exists():
        assert json.loads(path.read_text()) == report, 'Existing audit must replay identically'
    else:
        path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(passed=True, new_native_workflows=3, new_cached_references=3,
                         interrupted_input_sha256_identical=True), indent=2))


if __name__ == '__main__':
    main()
