"""Finite sequential runner, with explicit failure state; no scientific reruns."""
import os
import platform
import subprocess
from sensitivity import *


if __name__ == '__main__':
    state = dict(pid=os.getpid(), python=platform.python_version(), numpy=np.__version__,
        scientific_source_sha256=sha(P/'sensitivity.py'), protocol_sha256=sha(P/'PROTOCOL.md'),
        baseline_commit='f7f6c5754d19f5c345ae9ce3b2f368beb45cb1cb', completed=[], stage='running')
    for row in registry():
        tag = row['tag']
        if (P/f'result-{tag}.json').exists():
            state['completed'].append(tag); continue
        state['current'] = tag; save(P/'execution-state.json', state)
        try:
            with (P/f'run-{tag}.log').open('x') as log:
                subprocess.run([sys.executable, str(P/'sensitivity.py'), '--tag', tag],
                    stdout=log, stderr=subprocess.STDOUT, check=True)
        except Exception as error:
            state.update(stage='failed', error=str(error)); save(P/'execution-state.json', state)
            raise
        state['completed'].append(tag); save(P/'execution-state.json', state)
    state.update(stage='ready_for_interpretation', current=None); save(P/'execution-state.json', state)
    print('All nine views complete; response propagation, independent verification and interpretation remain.')
