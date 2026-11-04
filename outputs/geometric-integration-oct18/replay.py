"""Read-only replay of every Source Data 71 case through the packaged routines."""
import ast
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal.geometric_uncertainty import geometric_ratio_set
from codedcal.uncertainty import fieller_set

RESEARCH = ROOT/'outputs/temporal-variants-oct17'
SOURCE = ROOT/'outputs/temporal-validation-oct09'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):
            h.update(block)
    return h.hexdigest()


def main():
    start=time.monotonic()
    frozen=json.loads((RESEARCH/'freeze.json').read_text())
    for name,digest in frozen['files'].items():
        assert sha(ROOT/name)==digest,name
    sources=json.loads((RESEARCH/'source-hashes.json').read_text())
    counts=dict(arithmetic=0,geometric=0)
    covered={v:[0,0] for v in counts}
    reference_artifacts={}
    worlds=0;fields=0
    for file in sorted((RESEARCH/'worlds').glob('*.jsonl.gz')):
        marker=json.loads(file.with_suffix('.json').read_text())
        assert sha(file)==marker['output_sha256']
        reference_artifacts[str(file.relative_to(ROOT))]=marker['output_sha256']
        label=file.name.removesuffix('.jsonl.gz')
        rep=None;arrays=None
        with gzip.open(file,'rt') as stream:
            for line in stream:
                row=json.loads(line)
                if rep!=row['rep']:
                    rep=row['rep'];path=SOURCE/'worlds'/label/f'field-{rep:03d}.npz'
                    assert sha(path)==sources[str(path.relative_to(SOURCE))],str(path)
                    with np.load(path) as data:
                        arrays={k:data[k] for k in ('coded','known_source')}
                    fields+=1
                data=(arrays['coded'][row['assignment']] if row['method']=='coded' else arrays['known_source'])
                groups=data[:row['group_count'],0 if row['q']==.1 else 1]
                if row['variant']=='geometric':
                    new=geometric_ratio_set(groups)
                    for k in ('ratio','usable','kind','intervals','positive_reference'):
                        assert new[k]==row[k],(label,rep,k,new[k],row[k])
                    for actual,expected in zip(new['windows'],row['windows']):
                        assert all(actual[k]==value for k,value in expected.items()),(label,rep,'window')
                else:
                    pooled=groups.mean(axis=1)
                    new=fieller_set(pooled[:,0],pooled[:,1])
                    assert (new['ratio'] if new['positive_reference'] else None)==row['ratio']
                    for k in ('usable','kind','intervals','positive_reference'):
                        assert new[k]==row[k],(label,rep,k)
                counts[row['variant']]+=1
                if row['world']=='noise_shift' and row['scope']=='observed' and row['method']=='coded':
                    covered[row['variant']][0]+=1
                    covered[row['variant']][1]+=bool(new['usable'] and row['covers'])
        worlds+=1
        if worlds%14==0:print('Replayed',worlds,'worlds;',sum(counts.values()),'cases',flush=True)
    assert worlds==126 and fields==3528 and counts==dict(arithmetic=42336,geometric=42336)
    assert covered==dict(arithmetic=[864,627],geometric=[864,837])
    # Historical fitting functions must remain literally the same AST bodies.
    historical=subprocess.check_output(['git','show','374e966:tool/src/codedcal/smlm.py'],cwd=ROOT,text=True)
    current=(ROOT/'tool/src/codedcal/smlm.py').read_text()
    names=('cubic','fit_ratio','components','_fit','split_summary')
    def functions(source):
        return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(source).body
                if isinstance(n,ast.FunctionDef) and n.name in names}
    assert functions(historical)==functions(current)
    result=dict(passed=True,worlds=worlds,retained_fields=fields,cases=counts,
                geometric_window_sets_verified=2*counts['geometric'],
                equality='Exact equality of points, full interval unions, support flags and both window sets',
                uniform_kernel_joint_usable_coverage=covered,
                historical_fit_functions_unchanged=list(names),
                source71_freeze_sha256=sha(RESEARCH/'freeze.json'),
                source71_frozen_dependencies_verified=True,
                source71_case_artifacts=reference_artifacts,
                package_source_hashes={str(p.relative_to(ROOT)):sha(p) for p in
                    [ROOT/'tool/src/codedcal/geometric_uncertainty.py',ROOT/'tool/src/codedcal/uncertainty.py']},
                physical_validation=False,new_observations_generated=False,elapsed_s=time.monotonic()-start)
    (HERE/'replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Exact package replay passed:',counts,flush=True)


if __name__=='__main__':main()
