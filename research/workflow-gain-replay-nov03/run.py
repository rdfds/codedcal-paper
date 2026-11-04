"""Single finite two-stage replay. Historical reproduction gates all new outcomes."""
import fcntl
import gc
import importlib.metadata
import os
import platform
import resource
import shutil
import signal
import subprocess
import time
import traceback
from unittest.mock import patch
from replay_common import *
from fixed_decoder import FixedDecoder,HistoricalGate,catalogs_reproduced,CosineOptimizedPixelDecoder


class Stop(BaseException):pass


def main():
    f=load(H/'freeze.json')
    assert subprocess.check_output(['git','show','HEAD:'+str((H/'freeze.json').relative_to(R))],cwd=R)==(H/'freeze.json').read_bytes(),'Commit freeze first'
    assert platform.python_version()==f['python']
    for name,value in f['environment'].items():assert importlib.metadata.version(name)==value,name
    for name,digest in f['code_sha256'].items():assert sha(R/name)==digest,name
    for name,digest in f['serval_sha256'].items():assert sha(R.parent/'serval-upstream.nosync'/name)==digest,name
    assert sha(R/f['original_freeze_path'])==f['original_freeze_sha256']
    for name,pin in f['preceding_seals'].items():assert sha(R/name)==pin['sha256']
    assert sha(R/f['input']['path'])==f['input']['sha256']
    assert shutil.disk_usage(R).free>=f['minimum_start_free_bytes'],'Insufficient start disk'
    out=H/'run';out.mkdir(exist_ok=True)
    locks=[]
    for path in [P/'run/collector.lock',R/'research/workflow-failure-audit-nov03/native-run/collector.lock',out/'collector.lock']:
        lock=path.open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);locks.append(lock)
    assert not (out/'state.json').exists(),'Never restart or overwrite this replay'
    begin=time.monotonic();wall=time.time()
    state=dict(pid=os.getpid(),command=str(Path(__file__).resolve()),stage='starting',started_epoch=wall,
        freeze_sha256=sha(H/'freeze.json'),decoder_attempts=0,decoder_predictions_started=0,complete_predictions=0,
        fixed_codebook_index_builds=0,registration_optimizations=0,gain_optimizations=0,exact_cache_hits=0)
    def save():
        state.update(elapsed_monotonic_seconds=time.monotonic()-begin,elapsed_wall_seconds=time.time()-wall,
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,free_bytes=shutil.disk_usage(R).free)
        write_json(out/'state.json',state)
    def check():
        if time.monotonic()-begin>=f['maximum_wall_seconds']:raise Stop('Finite time limit')
        if shutil.disk_usage(R).free<f['minimum_running_free_bytes']:raise Stop('Disk reserve')
        if sum(p.stat().st_size for p in out.rglob('*') if p.is_file())>f['maximum_retained_bytes']:raise Stop('Output limit')
    def interrupted(signum,frame):raise Stop('Supervisor/time signal '+str(signum))
    for sig in [signal.SIGTERM,signal.SIGALRM]:signal.signal(sig,interrupted)
    signal.alarm(f['maximum_wall_seconds']);save()
    try:
        case=pad_for_challenges(ImageCase.load(R/f['input']['path']))
        engine=FixedDecoder(case,f['adapter_config']['source_native_scales']);state['fixed_codebook_index_builds']=1
        gate=HistoricalGate([e['id'] for e in f['cases']])
        native={p:load(P/'run/trials'/p/'native.json') for p in f['policies']}
        # An accidental optimizer call is an execution failure, never a scored molecular failure.
        with patch.object(CosineOptimizedPixelDecoder,'fit',side_effect=AssertionError('No optimization authorized')):
            for mode in ['historical','anchored']:
                if mode=='anchored':
                    gate.require_complete()
                    write_json(out/'historical-gate.json',dict(passed=True,cases=26,all_original_assessments_exact=True,
                        anchored_predictions_started=0,coordinate_tolerance=f['coordinate_engineering_tolerance']))
                for entry in f['cases']:
                    check();state.update(stage=mode,current_case=entry['id']);save()
                    assert sha(R/entry['source_path'])==entry['source_sha256']
                    for c in entry['matching_completed_calls']:assert sha(R/c['path'])==c['sha256']
                    raw=observed(case,entry['challenge'])
                    assert cache_key(raw,f['adapter_config'],entry['policy'])==entry['input_fingerprint']
                    old=retained(entry);maps=old['corrections'];gains=old['diagnostics']['fitted_source_gains']
                    path=out/mode/(entry['id']+'.json');assert not path.exists()
                    row=dict(stage='running',mode=mode,case_id=entry['id'],policy=entry['policy'],
                        challenge=entry['challenge'],input_fingerprint=entry['input_fingerprint'],
                        original_record_sha256=entry['source_sha256'],started_epoch=time.time())
                    state['decoder_attempts']+=1;save();write_json(path,row);start=time.perf_counter()
                    try:
                        result=engine.run(raw,maps,gains,mode)
                        assessment=score(result,native[entry['policy']],case,entry['challenge'],native[entry['policy']]['corrections'])
                        row.update(stage='completed',**result,assessment=assessment,
                            quality=blank_bound(Catalog.from_dict(result['catalog']),case),seconds=time.perf_counter()-start)
                        if mode=='historical':
                            matched=catalogs_reproduced(result['catalog'],old['catalog'],f['coordinate_engineering_tolerance']) and result['coverage']==old['coverage']
                            if entry['challenge'] is not None:
                                previous=load(R/entry['source_path'])['assessment'];matched=matched and assessment==previous
                            row['historical_reproduced']=bool(matched);gate.record(entry['id'],matched)
                            if not matched:raise ValueError('Historical reproduction failed: '+entry['id'])
                        write_json(path,row);state['complete_predictions']+=1;save()
                    except BaseException as exc:
                        row.update(stage='failed',error=type(exc).__name__+': '+str(exc),traceback=traceback.format_exc(),seconds=time.perf_counter()-start)
                        write_json(path,row);raise
                    finally:
                        state['decoder_predictions_started']=engine.predictions;save()
                    print(mode,entry['id'],assessment.get('failed'),assessment.get('identity_loss'),flush=True)
                    del raw,result,row;gc.collect()
        assert state['complete_predictions']==state['decoder_predictions_started']==52
        state['stage']='completed';save()
    except Stop as exc:
        state.update(stage='stopped_at_limit',error=str(exc));save()
    except BaseException as exc:
        state.update(stage='stopped_with_error',error=type(exc).__name__+': '+str(exc),traceback=traceback.format_exc());save();raise
    finally:signal.alarm(0)


if __name__=='__main__':main()
