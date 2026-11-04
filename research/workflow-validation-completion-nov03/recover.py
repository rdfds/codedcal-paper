"""One finite continuation: only the three missing frozen group-search probes."""
from pathlib import Path
import copy
import fcntl
import importlib.metadata
import json
import os
import shutil
import signal
import subprocess
import sys
import time
import traceback
import numpy as np

H = Path(__file__).resolve().parent
R = H.parents[1]
P = R / 'research/workflow-validation-controller-nov02'
sys.path[:0] = [str(P), str(R/'tool/src'), str(R/'research/workflow-memory-repair-nov02')]
from common import sha, cache_key
from codedcal.workflow_validation import ImageCase, Catalog, WorkflowResult, pad_for_challenges, run_challenge, compare, write_json


class Stop(BaseException):
    pass


def missing_pairs(freeze, run):
    return [(policy, ch) for policy in freeze['policies']
            if json.loads((run/'qualification'/(policy+'.json')).read_text())['primary_native_eligible']
            for ch in freeze['challenges'] if ch['id'] not in freeze['fractional']
            and not (run/'trials'/policy/(ch['id']+'.json')).exists()]


def restore_cache(records):
    """Cache completed observed outputs only; timing never counts as new work."""
    cache = {}
    for row in records:
        if row['stage'] != 'completed' or row['policy'] != 'group_search':
            continue
        key = row['input_fingerprint']
        value = WorkflowResult(Catalog.from_dict(row['catalog']), np.asarray(row['corrections']),
                               row['coverage'], copy.deepcopy(row['diagnostics']))
        if key in cache:
            old = cache[key]
            if old.catalog.as_dict() != value.catalog.as_dict() or not np.array_equal(old.corrections, value.corrections) or old.coverage != value.coverage:
                raise ValueError('Conflicting completed cache payload at fingerprint '+key)
        cache[key] = value
    return cache


def validate_dispatch(freeze, new_index, challenge, key, cached):
    if not 0 <= new_index < freeze['maximum_new_dispatches']:
        raise Stop('Finite six-dispatch recovery allowance exhausted')
    expected = freeze['missing_pairs'][new_index // 2]
    if expected['policy'] != 'group_search' or challenge != expected['challenge']['id']:
        raise Stop('Attempted an unfrozen workflow/challenge')
    if new_index == 0 and key != freeze['retry']['input_fingerprint']:
        raise Stop('Interrupted-call retry fingerprint changed')
    if new_index % 2 == 1 and not cached:
        raise Stop('Inverse reference is not in the verified exact result cache')


def main():
    freeze_path = H/'freeze.json'
    freeze = json.loads(freeze_path.read_text())
    committed = subprocess.check_output(['git', 'show', 'HEAD:'+str(freeze_path.relative_to(R))], cwd=R)
    assert committed == freeze_path.read_bytes(), 'Commit the recovery freeze first'
    assert sha(P/'freeze.json') == freeze['original_freeze_sha256']
    for name, digest in freeze['code_sha256'].items():
        assert sha(R/name) == digest, name
    for name, digest in freeze['recovery_code_sha256'].items():
        assert sha(R/name) == digest, name
    for key, dirname in [('bardensr_sha256', 'bardensr-upstream.nosync'), ('serval_sha256', 'serval-upstream.nosync')]:
        for name, digest in freeze[key].items():
            assert sha(R.parent/dirname/name) == digest, name
    for name, expected in freeze['environment'].items():
        assert importlib.metadata.version(name) == expected, name
    assert shutil.disk_usage(R).free >= freeze['minimum_start_free_bytes']
    out = H/'run'
    oldlock = (P/'run/collector.lock').open('r')
    fcntl.flock(oldlock, fcntl.LOCK_EX|fcntl.LOCK_NB)
    lock = (out/'collector.lock').open('a+')
    fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert not (out/'state.json').exists(), 'Never overwrite/restart this recovery'
    for row in json.loads((H/'prefix-manifest.json').read_text())['files']:
        assert sha(R/row['path']) == row['sha256'] == sha(R/row['source'])
    pairs = missing_pairs(freeze, out)
    assert [{'policy': p, 'challenge': ch} for p, ch in pairs] == freeze['missing_pairs']
    assert sha(R/freeze['input']['path']) == freeze['input']['sha256']
    case = pad_for_challenges(ImageCase.load(R/freeze['input']['path']))
    from memory_adapter import isolated
    base = isolated(freeze['adapter_config'], out/'children')
    prior_calls = [json.loads(p.read_text()) for p in sorted((out/'calls').glob('*.json'))]
    assert len(prior_calls) == 58
    base.cache.update(restore_cache(prior_calls))
    native_row = json.loads((out/'trials/group_search/native.json').read_text())
    native_key = cache_key(case.for_adapter(), freeze['adapter_config'], 'group_search')
    native = base.cache[native_key]
    assert native.catalog.as_dict() == native_row['catalog']
    assert np.array_equal(native.corrections, native_row['corrections'])
    write_json(out/'cache-import.json', dict(completed_group_records=sum(c['stage']=='completed' and c['policy']=='group_search' for c in prior_calls),
        unique_verified_results=len(base.cache), native_fingerprint=native_key, new_actual_workflows=0))
    begin_wall = time.time()
    begin_mono = time.monotonic()
    deadline = min(freeze['absolute_deadline_epoch'], begin_wall+freeze['wall_seconds_limit'])
    state = dict(pid=os.getpid(), command=str(Path(__file__).resolve()), stage='starting',
        started_epoch=begin_wall, deadline_epoch=deadline, freeze_sha256=sha(freeze_path),
        logical_call_starts=58, complete_calls=57, completed_pairs=21, inherited_call_count=58,
        completed_new_pairs=0, total_reserved=freeze['total_reserved'])

    def save():
        state.update(elapsed_wall_seconds=time.time()-begin_wall, elapsed_monotonic_seconds=time.monotonic()-begin_mono)
        write_json(out/'state.json', state)

    def check():
        if time.time() >= deadline or time.monotonic()-begin_mono >= freeze['wall_seconds_limit']:
            raise Stop('Bounded recovery deadline reached')
        if shutil.disk_usage(R).free < freeze['minimum_running_free_bytes']:
            raise Stop('Disk reserve reached')
        if sum(p.stat().st_size for p in out.rglob('*') if p.is_file() and p.name != 'input.npz') > freeze['maximum_output_bytes']:
            raise Stop('Finite output budget reached')

    def interrupted(signum, frame):
        raise Stop('Supervisor/time signal '+str(signum))

    for sig in (signal.SIGTERM, signal.SIGALRM):
        signal.signal(sig, interrupted)
    signal.alarm(max(1, int(deadline-time.time())))
    save()
    try:
        class Recorded:
            def context(self, c):
                return base.context(c)

            def run(self, c, policy):
                check()
                assert policy == 'group_search'
                index = state['logical_call_starts']
                key = cache_key(c, freeze['adapter_config'], policy)
                validate_dispatch(freeze, index-58, state['current_challenge'], key, key in base.cache)
                state['logical_call_starts'] += 1
                state.update(active_call=index, active_policy=policy)
                save()
                path = out/'calls'/f'{index:04}.json'
                assert not path.exists()
                start = time.perf_counter()
                entry = dict(stage='running', call=index, policy=policy, phase='frozen_probe_completion',
                    challenge=state['current_challenge'], started_epoch=time.time(), input_fingerprint=key)
                if index == 58:
                    entry['retry_of'] = freeze['retry']
                write_json(path, entry)
                try:
                    result = base.run(c, policy)
                    entry.update(stage='completed', catalog=result.catalog.as_dict(), corrections=result.corrections,
                        coverage=result.coverage, diagnostics=result.diagnostics, wall_seconds=time.perf_counter()-start)
                    write_json(path, entry)
                    state['complete_calls'] += 1
                    save()
                    return result
                except BaseException as exc:
                    entry.update(stage='failed', error=type(exc).__name__+': '+str(exc), traceback=traceback.format_exc(),
                                 wall_seconds=time.perf_counter()-start)
                    write_json(path, entry)
                    save()
                    raise

        adapter = Recorded()
        for policy, challenge in pairs:
            state.update(stage='comparative_validation', current_policy=policy, current_challenge=challenge['id'])
            save()
            check()
            path = out/'trials'/policy/(challenge['id']+'.json')
            assert not path.exists()
            row = run_challenge(case, adapter, policy, challenge, minimum_reference=freeze['reference_minimum'])
            reference = row['runs'][1]
            if reference['status'] == 'completed':
                row['oracle_vs_native'] = compare(native.catalog, Catalog.from_dict(reference['catalog']), case.blank_classes)
                exact = reference['catalog'] == native.catalog.as_dict() and np.array_equal(reference['corrections'], native.corrections)
                row['integer_oracle_native_exact'] = bool(exact)
                if not exact:
                    row['raw_assessment_before_control'] = row['assessment']
                    row['assessment'] = dict(status='unsupported', failed=None, reason='integer_oracle_native_mismatch')
            write_json(path, row)
            state['completed_pairs'] += 1
            state['completed_new_pairs'] += 1
            save()
            print(policy, challenge['label'], row['assessment'], flush=True)
        assert state['logical_call_starts'] == 64 and state['completed_pairs'] == 24
        state['stage'] = 'completed'
        save()
    except Stop as exc:
        state.update(stage='stopped_at_limit', error=str(exc))
        save()
        print('STOP', str(exc), flush=True)
    except BaseException as exc:
        state.update(stage='stopped_with_error', error=type(exc).__name__+': '+str(exc), traceback=traceback.format_exc())
        save()
        raise
    finally:
        signal.alarm(0)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
