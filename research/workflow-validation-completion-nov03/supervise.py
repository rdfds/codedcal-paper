"""Use the existing finite process-group supervisor with a new explicit cap."""
from pathlib import Path
import json
import sys
import time

H = Path(__file__).resolve().parent
R = H.parents[1]
sys.path.insert(0, str(R/'research/workflow-memory-repair-nov02'))
from supervise_finite import main

if __name__ == '__main__':
    f = json.loads((H/'freeze.json').read_text())
    remaining = int(f['absolute_deadline_epoch']-time.time())-2
    assert remaining > 0
    sys.argv = [sys.argv[0], '--worker', str(H/'recover.py'), '--state', str(H/'supervisor-state.json'),
                '--seconds', str(min(f['wall_seconds_limit'], remaining))]
    raise SystemExit(main())
