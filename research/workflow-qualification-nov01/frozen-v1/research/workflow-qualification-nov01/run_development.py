"""One finite, exclusive, committed-freeze development benchmark."""
from pathlib import Path
import os,sys,json,hashlib,time,signal,fcntl,traceback,subprocess
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path.insert(0,str(H))
from native_adapters import make_adapter
from codedcal.workflow_validation import ImageCase,Catalog,WorkflowResult,pad_for_challenges,run_challenge,write_json,compare
F=H/'development-freeze.json';freeze=json.loads(F.read_text());OUT=H/'development';OUT.mkdir(exist_ok=True)
lock=(OUT/'collector.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (OUT/'state.json').exists(),'No silent restart/overwrite; diagnose retained evidence first'
for p,d in freeze['code_sha256'].items():assert hashlib.sha256((R/p).read_bytes()).hexdigest()==d,p
for p,d in freeze['native_code_sha256'].items():assert hashlib.sha256((R.parent/'bardensr-upstream.nosync'/p).read_bytes()).hexdigest()==d,p
assert subprocess.check_output(['git','show','HEAD:'+str(F.relative_to(R))],cwd=R)==F.read_bytes(),'Commit freeze before outcomes'
start=time.monotonic();state=dict(pid=os.getpid(),command=str(Path(__file__).resolve()),stage='running',started_epoch=time.time(),completed_pairs=0,logical_charged=0,freeze_sha256=hashlib.sha256(F.read_bytes()).hexdigest())
class Deadline(BaseException):pass
signal.signal(signal.SIGALRM,lambda *args: (_ for _ in ()).throw(Deadline()));signal.alarm(freeze['worker_seconds_limit'])
def save():state['elapsed_seconds']=time.monotonic()-start;write_json(OUT/'state.json',state)
save()
try:
    # Source/backend/control exclusion is fixed before outcomes, never chosen by a method win.
    for e in freeze['cases']:
        case=pad_for_challenges(ImageCase.load(e['path']));assert hashlib.sha256(Path(e['path']).read_bytes()).hexdigest()==e['sha256']
        geometry_cache={}
        for decoder in e['decoders']:
            if not e['qualified'][decoder]:continue
            adapter=make_adapter({'decoder':decoder});adapter.geometry_cache=geometry_cache
            for policy in freeze['policies']:
                state.update(current_case=e['name'],current_decoder=decoder,current_policy=policy,current_challenge='native');save()
                dest=OUT/e['name']/decoder/policy;dest.mkdir(parents=True,exist_ok=True)
                state['logical_charged']+=47;assert state['logical_charged']<=freeze['logical_limit'];save();begin=time.perf_counter()
                try:
                    result=adapter.run(case.for_adapter(),policy)
                    native=dict(status='completed',catalog=result.catalog.as_dict(),corrections=result.corrections,coverage=result.coverage,diagnostics=result.diagnostics,wall_seconds=time.perf_counter()-begin)
                except Exception as x:native=dict(status='failed',error=type(x).__name__+': '+str(x),traceback=traceback.format_exc(),wall_seconds=time.perf_counter()-begin)
                write_json(dest/'native.json',native)
                for ch in e['plan']['challenges']:
                    state['current_challenge']=ch['id'];save()
                    path=dest/(ch['id']+'.json');assert not path.exists()
                    row=run_challenge(case,adapter,policy,ch)
                    # Native consistency of the composed reference is retained separately.
                    b=row['runs'][1]
                    if native['status']=='completed' and b['status']=='completed':row['oracle_vs_native']=compare(Catalog.from_dict(native['catalog']),Catalog.from_dict(b['catalog']),case.blank_classes)
                    write_json(path,row);state['completed_pairs']+=1;save()
                    print(e['name'],decoder,policy,ch['role'],row['assessment']['status'],row['assessment'].get('identity_loss'),flush=True)
                # Charge each strategy separately despite exact cross-strategy reuse.
                save()
        geometry_cache.clear()
    state['stage']='completed';save()
except Deadline:
    state.update(stage='time_limit',error='Finite three-hour deadline reached; current incomplete invocation and all completed evidence retained');save()
except BaseException as x:
    state.update(stage='stopped_with_error',error=type(x).__name__+': '+str(x),traceback=traceback.format_exc());save();raise
finally:signal.alarm(0)
