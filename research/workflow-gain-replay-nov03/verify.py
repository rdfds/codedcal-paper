"""Replay the recorded adjudication and provenance without decoding images."""
import argparse
import fcntl
import subprocess
from replay_common import *
from codedcal.workflow_gain import anchor_relative_gains


def main(seal=False):
    f=load(H/'freeze.json');out=H/'run'
    lock=(out/'collector.lock').open('r');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    state=load(out/'state.json');assert state['stage'] in {'completed','stopped_at_limit','stopped_with_error'}
    for path,digest in f['code_sha256'].items():assert sha(R/path)==digest,path
    for path,pin in f['preceding_seals'].items():
        assert sha(R/path)==pin['sha256'];assert verify_manifest(R/path)==pin['files_verified']
    assert sha(R/f['original_freeze_path'])==f['original_freeze_sha256']
    assert sha(R/f['input']['path'])==f['input']['sha256']
    for path,digest in f['serval_sha256'].items():assert sha(R.parent/'serval-upstream.nosync'/path)==digest,path
    # Metadata-only arrays are sufficient for re-scoring; no image array is read.
    from types import SimpleNamespace
    with np.load(R/f['input']['path'],allow_pickle=False) as a:
        case=SimpleNamespace(class_ids=a['class_ids'],blank_classes=a['blank_classes'])
    from fixed_decoder import catalogs_reproduced
    completed=0;baseline=0;started=0
    for e in f['cases']:
        assert sha(R/e['source_path'])==e['source_sha256']
        for p in e['matching_completed_calls']:assert sha(R/p['path'])==p['sha256']
        old=retained(e);native=load(P/'run/trials'/e['policy']/'native.json')
        for mode in ['historical','anchored']:
            path=out/mode/(e['id']+'.json')
            if not path.exists():continue
            r=load(path);started+=1
            if r['stage']!='completed':continue
            completed+=1
            assert r['input_fingerprint']==e['input_fingerprint'] and r['corrections']==old['corrections']
            gains=np.asarray(old['diagnostics']['fitted_source_gains'])
            if mode=='anchored':
                assert load(out/'historical-gate.json')['passed']
                reference=np.asarray(f['adapter_config']['source_native_scales']);reference/=reference.mean()
                expected,gauge=anchor_relative_gains(gains,reference)
                assert np.array_equal(expected,r['gains_used']) and gauge==r['gain_gauge']
            else:
                assert np.array_equal(gains,r['gains_used']) and r['gain_gauge'] is None
                assert catalogs_reproduced(r['catalog'],old['catalog']) and r['coverage']==old['coverage'] and r['historical_reproduced']
                baseline+=1
            assert r['assessment']==score(r,native,case,e['challenge'],native['corrections'])
            assert r['quality']==blank_bound(Catalog.from_dict(r['catalog']),case)
            if mode=='historical' and e['challenge'] is not None:assert r['assessment']==load(R/e['source_path'])['assessment']
            assert r['decoder_only_executions']==1 and r['calibration_optimizations']==r['gain_optimizations']==0
    assert completed==state['complete_predictions'] and started==state['decoder_attempts']
    if state['stage']=='completed':assert completed==52 and baseline==26
    # Unlock before the separate read-only analyzer takes the same exclusive lock.
    fcntl.flock(lock,fcntl.LOCK_UN)
    subprocess.run([sys.executable,str(H/'analyze.py')],check=True)
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    manifest=H/'completion-manifest.json'
    if seal:
        assert not manifest.exists(),'Never overwrite sealed evidence'
        paths=[p for p in H.rglob('*') if p.is_file() and '__pycache__' not in p.parts
            and p.name not in {'completion-manifest.json','verification.json'}]
        write_json(manifest,dict(files=[dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(paths)]))
    n=verify_manifest(manifest) if manifest.exists() else None
    result=dict(passed=True,completed_decoder_predictions=completed,historical_cases_reproduced=baseline,
        all_24_original_assessments_unchanged=baseline==26,original_references_and_parameters_verified=True,
        anchored_gain_only_intervention_verified=True,preceding_seals_unchanged=f['preceding_seals'],
        sealed_files=n,new_predictions_during_verification=0,physical_validation_added=False)
    write_json(H/'verification.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--seal',action='store_true');main(parser.parse_args().seal)
