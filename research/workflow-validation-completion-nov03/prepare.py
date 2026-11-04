"""Verify sealed evidence and freeze a checkpoint-only bounded recovery."""
from pathlib import Path
import copy
import datetime
import hashlib
import json
import shutil
import subprocess
import sys
import time

H = Path(__file__).resolve().parent
R = H.parents[1]
P = R/'research/workflow-validation-controller-nov02'
sys.path.insert(0, str(H))
from recover import sha, missing_pairs, restore_cache


def verify_manifest(path):
    data = json.loads(path.read_text())
    for row in data['files']:
        p = R/row['path']
        assert p.is_file() and not p.is_symlink() and p.stat().st_size == row['bytes'], row['path']
        assert sha(p) == row['sha256'], row['path']
    return len(data['files'])


def main():
    assert not (H/'freeze.json').exists()
    f = json.loads((P/'freeze.json').read_text())
    state = json.loads((P/'run/state.json').read_text())
    assert state['stage'] == 'stopped_at_limit' and state['completed_pairs'] == 21
    assert state['logical_call_starts'] == 58 and state['complete_calls'] == 57
    counts = {}
    for name in ['workflow-validation-controller-nov02/closure-manifest.json',
                 'workflow-validation-nov02/completion-manifest.json',
                 'workflow-memory-repair-nov02/repair-manifest.json',
                 'workflow-validation-memory-nov02/stopped-manifest.json',
                 'workflow-validation-checkpoint-nov02/stopped-manifest.json']:
        counts[name] = verify_manifest(R/'research'/name)
    for path, digest in f['code_sha256'].items():
        assert sha(R/path) == digest, path
    assert sha(R/f['input']['path']) == f['input']['sha256']
    missing = missing_pairs(f, P/'run')
    summary = json.loads((P/'summary.json').read_text())
    assert [p+'/'+c['id'] for p, c in missing] == summary['missing']
    assert len(missing) == 3 and all(p == 'group_search' and ch['role'] == 'probe' for p, ch in missing)
    assert summary['adequate_cases'] == 0 and summary['failed_cases'] == 2
    calls = [json.loads(p.read_text()) for p in sorted((P/'run/calls').glob('*.json'))]
    assert [c['call'] for c in calls] == list(range(58))
    partial = [c for c in calls if c['stage'] != 'completed']
    assert len(partial) == 1 and partial[0]['call'] == 57
    assert partial[0]['challenge'] == missing[0][1]['id']
    cache = restore_cache(calls)
    assert cache and partial[0]['input_fingerprint'] not in cache

    rows = []
    for directory in ['calls', 'qualification', 'trials']:
        for source in sorted((P/'run'/directory).rglob('*.json')):
            target = H/'run'/source.relative_to(P/'run')
            assert not target.exists()
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            digest = sha(source)
            assert sha(target) == digest
            rows.append(dict(path=str(target.relative_to(R)), source=str(source.relative_to(R)),
                             bytes=target.stat().st_size, sha256=digest))
    prefix = dict(files=rows, inherited_dispatches=58, inherited_complete_calls=57,
        inherited_pairs=21, inherited_actual_workflows=26, inherited_exact_cache_hits=31,
        inherited_interrupted_calibrations=1, new_fits=0,
        previous_sealed_files_verified=counts, original_closure_manifest_sha256=sha(P/'closure-manifest.json'))
    (H/'prefix-manifest.json').write_text(json.dumps(prefix, indent=2)+'\n')

    # Decision calculations stay byte-identical. Change only three report strings
    # whose original nine-call/historical labels no longer describe this prefix.
    text = (P/'analyze.py').read_text()
    replacements = {
        '# Repaired-selector prospective source check': '# Completion of the frozen public-source qualification check',
        'including the nine inherited checkpoints': 'including all 58 inherited dispatches and the preserved interruption',
        'Native controls imported from the old record are not newly measured runs': 'All inherited checkpoints and imported native controls are not newly measured runs',
    }
    for a, b in replacements.items():
        assert text.count(a) == 1
        text = text.replace(a, b)
    (H/'analyze.py').write_text(text)
    f = copy.deepcopy(f)
    f.update(inherited_call_count=58, total_reserved=1151, original_reservation_preserved=1150,
        additional_interrupted_retry_reserved=1, maximum_new_dispatches=6,
        missing_pairs=[dict(policy=p, challenge=ch) for p, ch in missing],
        original_freeze_sha256=sha(P/'freeze.json'), original_closure_manifest_sha256=sha(P/'closure-manifest.json'),
        continuation_of=str(P.relative_to(R)), prior_inherited_native_count=9,
        original_deadline_epoch=f['continuation_deadline_epoch'],
        absolute_deadline_epoch=time.time()+5400, wall_seconds_limit=3600,
        minimum_start_free_bytes=4*1024**3,
        user_authorization='After cache cleanup, user: ok so finish it',
        source_comparative_outcomes_examined=True, comparative_outcomes_examined=True,
        prior_partial_result=dict(adequate=0, failed=2, conventional_flags=2, weighted_flags=1, weighted_unresolved=1),
        retry=dict(path=str((P/'run/calls/0057.json').relative_to(R)), sha256=sha(P/'run/calls/0057.json'),
                   input_fingerprint=partial[0]['input_fingerprint'], original_call=57, maximum_retries=1),
        analysis_source_sha256=sha(P/'analyze.py'), analysis_text_only_replacements=replacements,
        preceding_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=R, text=True).strip(),
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    own = list(H.glob('*.py'))+[H/'PROTOCOL.md', H/'prefix-manifest.json']
    f['recovery_code_sha256'] = {str(p.relative_to(R)): sha(p) for p in sorted(own)}
    assert f['total_reserved'] == 1151 <= f['maximum_total'] == 1600
    (H/'freeze.json').write_text(json.dumps(f, indent=2)+'\n')
    print(json.dumps(dict(verified_prior_manifests=counts, prefix_files=len(rows), inherited_pairs=21,
        missing_probe_pairs=3, maximum_new_dispatches=6, total_reserved=1151,
        new_actual_workflows=0, verified_cache_keys=len(cache)), indent=2))


if __name__ == '__main__':
    main()
