"""One finite subprocess per actual workflow; exact result caches stay in parent.

This changes process lifetime, not calibration/decoder parameters. Local input
copies are bounded, hashed and removed only after proving they are unchanged
redundant serializations. Child results, logs and resource receipts are retained.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import importlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
import traceback
import numpy as np
from .workflow_validation import ImageCase, Catalog, WorkflowResult, write_json
from .workflow_bardensr import BardensrAdapter


def file_sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def cache_key(case, config, policy):
    h = hashlib.sha256()
    for a in [case.images, case.support, case.codebook, case.rounds, case.colors, case.class_ids, case.blank_classes]:
        a = np.ascontiguousarray(a)
        h.update(str((a.dtype.str, a.shape)).encode())
        h.update(memoryview(a).cast('B'))
    h.update(json.dumps({k: v for k, v in config.items() if k != 'decoder'}, sort_keys=True).encode())
    h.update(json.dumps(case.metadata, sort_keys=True).encode())
    h.update(policy.encode())
    return h.hexdigest()


class IsolatedWorkflowAdapter:
    def __init__(self, config, backend, workdir, maximum_input_bytes=64*1024**2):
        self.config = dict(config)
        self.backend = backend
        self.workdir = Path(workdir)
        self.maximum_input_bytes = maximum_input_bytes
        self.cache = {}
        self.geometry_cache = {}
        self.sequence = 0

    context = BardensrAdapter.context

    def run(self, case, policy):
        start = time.perf_counter()
        key = cache_key(case, self.config, policy)
        if key in self.cache:
            result = copy.deepcopy(self.cache[key])
            for name in ['actual_calibration_executions', 'actual_decoder_executions', 'actual_workflow_executions',
                         'child_process_executions', 'child_peak_rss_bytes']:
                result.diagnostics[name] = 0
            for name in ['calibration_seconds', 'decoder_seconds', 'child_workflow_seconds']:
                result.diagnostics[name] = 0.
            result.diagnostics.update(exact_result_cache_hit=True, wall_seconds=time.perf_counter()-start)
            return result
        arrays = {name: getattr(case, name) for name in
                  ['images', 'support', 'codebook', 'rounds', 'colors', 'class_ids', 'blank_classes']}
        if sum(a.nbytes for a in arrays.values()) + 65536 > self.maximum_input_bytes:
            raise ValueError('Finite local input serialization budget exceeded')
        self.workdir.mkdir(parents=True, exist_ok=True)
        job = self.workdir / f'{self.sequence:04}'
        self.sequence += 1
        job.mkdir(exist_ok=False)
        source = job / 'input.npz'
        np.savez(source, **arrays, name=case.name, preparation=case.preparation,
                 metadata_json=json.dumps(case.for_adapter().metadata))
        assert source.stat().st_size <= self.maximum_input_bytes
        config = job / 'config.json'
        write_json(config, self.config)
        output = job / 'result.json'
        receipt = dict(stage='prepared', policy=policy, input_fingerprint=key,
                       input_sha256=file_sha(source), local_input_bytes=source.stat().st_size,
                       backend=self.backend, started_epoch=time.time(), input_removed=False)
        path = job / 'receipt.json'
        write_json(path, receipt)
        process = None
        try:
            command = [sys.executable, '-m', 'codedcal.workflow_isolated', '--backend', self.backend,
                       '--input', str(source), '--config', str(config), '--policy', policy, '--output', str(output)]
            with (job / 'worker.log').open('x') as log:
                inherited = [p for p in os.environ.get('PYTHONPATH', '').split(os.pathsep) if p]
                paths = os.pathsep.join(dict.fromkeys([str(p) for p in sys.path if p] + inherited))
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                    env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONPATH=paths))
                receipt.update(stage='running', pid=process.pid, command=command)
                write_json(path, receipt)
                code = process.wait()
            payload = json.loads(output.read_text()) if output.exists() else {}
            if code or payload.get('stage') != 'completed':
                raise RuntimeError('Isolated native workflow failed: ' + str(payload.get('error', code)))
            result = WorkflowResult(Catalog.from_dict(payload['catalog']), np.asarray(payload['corrections']),
                                    payload['coverage'], payload['diagnostics'])
            result.diagnostics.update(child_workflow_seconds=result.diagnostics['wall_seconds'],
                child_process_executions=1, child_peak_rss_bytes=payload['peak_rss_bytes'],
                isolated_native_workflow=True, wall_seconds=time.perf_counter()-start)
            self.cache[key] = copy.deepcopy(result)
            receipt.update(stage='completed', exit_code=code, result_sha256=file_sha(output),
                           child_peak_rss_bytes=payload['peak_rss_bytes'])
            return result
        except BaseException as error:
            if process is not None and process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            receipt.update(stage='failed_or_interrupted', error=type(error).__name__ + ': ' + str(error),
                           traceback=traceback.format_exc(), exit_code=process.returncode if process else None)
            raise
        finally:
            # The source case + frozen challenge can regenerate these exact
            # bytes. Keep receipts/results/logs, not another full image copy.
            if source.exists():
                assert source.stat().st_size == receipt['local_input_bytes']
                assert file_sha(source) == receipt['input_sha256']
                source.unlink()
                receipt['input_removed'] = True
            receipt['wall_seconds'] = time.perf_counter()-start
            write_json(path, receipt)


def child_main():
    parser = argparse.ArgumentParser()
    for name in ['backend', 'input', 'config', 'policy', 'output']:
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    path = Path(args.output)
    try:
        module, factory = args.backend.split(':')
        adapter = getattr(importlib.import_module(module), factory)(json.loads(Path(args.config).read_text()))
        case = ImageCase.load(args.input).for_adapter()
        result = adapter.run(case, args.policy)
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if sys.platform != 'darwin':
            peak *= 1024
        write_json(path, dict(stage='completed', catalog=result.catalog.as_dict(), corrections=result.corrections,
                             coverage=result.coverage, diagnostics=result.diagnostics, peak_rss_bytes=peak))
        return 0
    except BaseException as error:
        write_json(path, dict(stage='failed', error=type(error).__name__ + ': ' + str(error),
                             traceback=traceback.format_exc()))
        raise


if __name__ == '__main__':
    raise SystemExit(child_main())
