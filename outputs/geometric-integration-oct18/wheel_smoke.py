"""Build/install offline and exercise the installed CLI outside the checkout."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
import tempfile

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def main():
    tmp=Path(tempfile.mkdtemp(prefix='codedcal-0.4-check-'))
    # Keep the venv executable path: resolving its symlink bypasses the venv.
    python=sys.executable
    log=[]
    def run(args, expected=0, env=None):
        result=subprocess.run([python]+args,cwd=tmp,env=env,text=True,capture_output=True)
        log.append('$ '+' '.join([python]+args)+'\n'+result.stdout+result.stderr)
        (HERE/'wheel-smoke.log').write_text('\n'.join(log))
        assert result.returncode==expected,(args,result.returncode,result.stdout,result.stderr)
        return result
    run(['-m','pip','wheel',str(ROOT/'tool'),'--no-index','--no-deps','--no-build-isolation','--wheel-dir',str(tmp/'wheels')])
    wheel=next((tmp/'wheels').glob('*.whl'))
    run(['-m','pip','install','--no-index','--no-deps','--target',str(tmp/'installed'),str(wheel)])
    env=dict(os.environ,PYTHONPATH=str(tmp/'installed'),PYTHONDONTWRITEBYTECODE='1')
    check=run(['-c','import codedcal; print(codedcal.__file__); print(codedcal.__version__)'],env=env)
    assert str(tmp/'installed') in check.stdout and '0.4.0' in check.stdout
    run(['-m','codedcal','--version'],env=env)
    events=np.array([[g,0,p,0] for g in range(8) for p in range(4)],float)
    np.savetxt(tmp/'events.csv',events,delimiter=',',header='group,code,panel,x',comments='')
    run(['-m','codedcal','smlm-budget',str(tmp/'events.csv'),'--uncertainty-variant','both',
         '--confidence-level','.9','--json',str(tmp/'supported.json')],env=env)
    np.savez(tmp/'weak.npz',z=np.zeros((8,4,4,1),complex),n=np.zeros((8,4,4)),q=[.1])
    run(['-m','codedcal','smlm-budget',str(tmp/'weak.npz'),'--uncertainty-variant','geometric',
         '--json',str(tmp/'weak.json')],expected=2,env=env)
    supported=json.loads((tmp/'supported.json').read_text())
    weak=json.loads((tmp/'weak.json').read_text())
    assert supported['reported_uncertainty_variants']==['arithmetic','geometric']
    assert supported['reported_uncertainty_usable']
    assert supported['geometric_uncertainty']['window_confidence_level']==.95
    assert weak['geometric_uncertainty']['per_frequency'][0]['kind']=='unbounded_positive'
    assert not weak['reported_uncertainty_usable']
    for name in ['supported.json','weak.json']:
        (HERE/('wheel-'+name)).write_bytes((tmp/name).read_bytes())
    record=dict(passed=True,version='0.4.0',wheel_name=wheel.name,
                wheel_sha256=hashlib.sha256(wheel.read_bytes()).hexdigest(),
                built_and_installed_without_network_or_dependencies=True,
                installed_to_isolated_target=True,executed_outside_source_checkout=True,
                supported_csv_both_exit=0,unresolved_npz_geometric_exit=2,
                json_retains_both_variants=True,physical_validation=False)
    (HERE/'wheel-smoke.json').write_text(json.dumps(record,indent=2)+'\n')
    shutil.copyfile(wheel,ROOT.parent.parent/'outputs'/wheel.name)
    print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':main()
