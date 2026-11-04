"""One sequential scientific worker; waits for verified existing inputs."""
import argparse
import fcntl
import json
import subprocess
import sys
import time
from common import CACHE,P,ROOT,cd20_registry,save

def main(skip):
    fields=[r['tag'] for r in cd20_registry()]
    # Available source projections allow both SDC views to be audited without
    # downloading another source. Simulations run only on the six CD20 fields.
    todo=['npc',*fields,'sdc_corrected','sdc_raw']
    state=P/'pipeline-state.json'
    for tag in todo:
        if tag in fields:
            row=next(r for r in cd20_registry() if r['tag']==tag)
            while True:
                try:manifest=json.loads((P/'input-manifest.json').read_text())
                except (FileNotFoundError,json.JSONDecodeError):manifest={}
                wanted=[row['prefix']+f'RESI/RESI rounds/RESI_R{r}_apicked.hdf5' for r in range(1,5)]
                if all(name in manifest for name in wanted):break
                save(state,dict(stage='waiting_verified_inputs',tag=tag))
                time.sleep(5)
        dest=P/f'characterization-{tag}.json'
        if not dest.exists():
            save(state,dict(stage='characterizing',tag=tag))
            subprocess.run([sys.executable,str(P/'reference_characterize.py'),'--tag',tag],check=True,cwd=ROOT)
        if tag in fields:
            dest=P/f'controls-{tag}.json'
            if not dest.exists() or not json.loads(dest.read_text())['complete']:
                save(state,dict(stage='simulating',tag=tag))
                subprocess.run([sys.executable,str(P/'simulate_reference.py'),'--tag',tag,'--draws','4000'],check=True,cwd=ROOT)
    save(state,dict(stage='ready_for_interpretation',distinct_acquisitions=8,views=9))

if __name__=='__main__':
    CACHE.mkdir(exist_ok=True)
    with (CACHE/'science.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        main(None)
