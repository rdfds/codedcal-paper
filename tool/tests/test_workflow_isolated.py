import json
import os
import subprocess
import time
from pathlib import Path
import numpy as np
import pytest
from codedcal.workflow_isolated import IsolatedWorkflowAdapter, cache_key
from codedcal.workflow_validation import ImageCase


def case():
    a = np.ones((4, 1, 8, 8), np.float32)
    return ImageCase(a, a.copy(), np.ones((4, 2)), np.array([0, 0, 1, 1]), np.array([0, 1, 0, 1]),
                     np.arange(2), np.array([1]), metadata=dict(challenge_truth='never send', context=[[0, 1], [0, 8], [0, 8]]))


def adapter(tmp_path, monkeypatch, maximum=64*1024**2):
    source = Path(__file__).resolve().parent
    monkeypatch.setenv('PYTHONPATH', os.environ.get('PYTHONPATH', '') + os.pathsep + str(source))
    return IsolatedWorkflowAdapter(dict(marker=str(tmp_path/'pids.txt')),
        'isolation_backend:make_adapter', tmp_path/'children', maximum)


def test_real_subprocess_cache_distinct_lifetime_and_redundant_input_cleanup(tmp_path, monkeypatch):
    a = adapter(tmp_path, monkeypatch)
    c = case().for_adapter()
    first = a.run(c, 'fixture')
    assert first.diagnostics['actual_workflow_executions'] == 1
    assert first.diagnostics['child_peak_rss_bytes'] > 0
    cached = a.run(c, 'fixture')
    assert cached.diagnostics['exact_result_cache_hit']
    assert cached.diagnostics['actual_workflow_executions'] == cached.diagnostics['child_peak_rss_bytes'] == 0
    assert cached.diagnostics['child_process_executions'] == 0
    changed = c.replace_images(c.images.copy(), c.support.copy())
    changed.images[0, 0, 0, 0] = 2
    a.run(changed, 'fixture')
    pids = (tmp_path/'pids.txt').read_text().splitlines()
    assert len(pids) == len(set(pids)) == 2 and str(os.getpid()) not in pids
    assert not list((tmp_path/'children').glob('*/input.npz'))
    for path in (tmp_path/'children').glob('*/receipt.json'):
        r = json.loads(path.read_text())
        assert r['stage'] == 'completed' and r['input_removed']
        assert len(r['input_sha256']) == len(r['result_sha256']) == 64


def test_child_failure_is_retained_and_never_cached(tmp_path, monkeypatch):
    a = adapter(tmp_path, monkeypatch)
    with pytest.raises(RuntimeError, match='Deliberate child failure'):
        a.run(case().for_adapter(), 'fail')
    receipt = json.loads((tmp_path/'children/0000/receipt.json').read_text())
    assert receipt['stage'] == 'failed_or_interrupted' and receipt['input_removed']
    assert (tmp_path/'children/0000/result.json').exists() and not a.cache


def test_oversize_input_does_not_launch_or_write(tmp_path, monkeypatch):
    a = adapter(tmp_path, monkeypatch, maximum=1)
    with pytest.raises(ValueError, match='serialization budget'):
        a.run(case().for_adapter(), 'fixture')
    assert not (tmp_path/'children').exists()


def test_interruption_terminates_and_reaps_child_before_propagating(tmp_path, monkeypatch):
    a = adapter(tmp_path, monkeypatch)
    original = subprocess.Popen
    class Interrupted(original):
        def wait(self, timeout=None):
            if timeout is None:
                time.sleep(.2)
                raise KeyboardInterrupt('Deliberate parent interruption')
            return super().wait(timeout=timeout)
    monkeypatch.setattr('codedcal.workflow_isolated.subprocess.Popen', Interrupted)
    with pytest.raises(KeyboardInterrupt):
        a.run(case().for_adapter(), 'wait')
    receipt = json.loads((tmp_path/'children/0000/receipt.json').read_text())
    assert receipt['stage'] == 'failed_or_interrupted' and receipt['input_removed']
    assert receipt['exit_code'] is not None
    with pytest.raises(ProcessLookupError):
        os.kill(receipt['pid'], 0)
