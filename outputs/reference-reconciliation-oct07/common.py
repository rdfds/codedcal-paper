"""Fixed exposed-cohort registry and utilities; no historical output mutation."""
from pathlib import Path
import hashlib
import json
import numpy as np

P = Path(__file__).resolve().parent
ROOT = P.parent.parent
CACHE = ROOT.parent / 'reference-reconciliation-oct07.nosync'

def input_path(row, r):
    """Resolve retained inputs in this checkout, preserving recorded provenance."""
    if row['tag'].startswith('sdc_'):
        return ROOT.parent/'heldout-nup96.nosync/independent-field'/f'round{r}_reduced.npz'
    base=ROOT.parent/'shaw-inputs.nosync/npc/extract_unrar' if row['tag']=='npc' else CACHE/'raw'
    return base/row['prefix']/f'RESI/RESI rounds/RESI_R{r}_apicked.hdf5'

def read(relative):
    return json.loads((ROOT / relative).read_text())

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while b := f.read(1 << 22):
            h.update(b)
    return h.hexdigest()

def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value

def save(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + '\n')

def cd20_registry():
    reps = read('outputs/replication-edf1112-sep15/field_results_edf1112.json')['reps']
    rows = []
    specs = [
        ('edf10rep1', 'EDF10/rep1/', 'edf10-replication-sep15',
         'field_results_edf10.json', 'pick_centers_edf10.npz', 'oracle_sources_edf10.npz'),
        ('edf10rep2', 'EDF10/rep2/', 'replication-cd20-sep15',
         'field_results_edf10rep2.json', 'pick_centers_edf10rep2.npz', 'oracle_sources_edf10rep2.npz'),
        ('fig4untreated', 'Figure4/untreated/', 'replication-cd20-sep15',
         'field_results_fig4untreated.json', 'pick_centers_fig4untreated.npz', 'oracle_sources_fig4untreated.npz'),
    ]
    for tag, prefix, folder, result, center, oracle in specs:
        base = ROOT / 'outputs' / folder
        old = read(str((base / result).relative_to(ROOT)))
        panels = old.get('panels_per_round', old.get('panels'))
        if isinstance(panels, list):
            panels = {str(r): panels for r in range(1, 5)}
        elif any(k.startswith('round_') for k in panels):
            panels = {k.replace('round_', ''): v for k, v in panels.items()}
        frames = old.get('frames_per_round', old.get('frames'))
        if isinstance(frames, int):
            frames = {str(r): frames for r in range(1, 5)}
        rows.append(dict(tag=tag, prefix=prefix, old=old, panels=panels,
                         frames={str(k): v for k, v in frames.items()},
                         center=str(base / center), oracle=str(base / oracle)))
    for rep, old in reps.items():
        base = ROOT / 'outputs/replication-edf1112-sep15'
        rows.append(dict(tag='rtx_' + rep, prefix=f'EDF11,EDF12/{rep}/',
                         old=old, panels=old['panels'], frames=old['frames'],
                         center=str(base / f'pick_centers_{rep}.npz'),
                         oracle=str(base / f'oracle_sources_{rep}.npz')))
    return rows

def historical_members():
    a = read('outputs/edf10-replication-sep15/fetch_manifest_edf10.json')['files']
    b = read('outputs/replication-cd20-sep15/fetch_manifest_cd20.json')['fields']
    c = read('outputs/replication-edf1112-sep15/fetch_manifest_edf1112.json')['reps']
    out = dict(a)
    for v in list(b.values()) + list(c.values()):
        out.update(v['files'])
    return out
