"""Offline closure; never fits workflows or contacts a data source."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

H=Path(__file__).resolve().parent
R=H.parents[1]
BASELINE='3789b2973740f6660846b951228095e6f9966496'


def digest(path,algorithm='sha256',prefix=b''):
    h=hashlib.new(algorithm);h.update(prefix)
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(4*1024**2),b''):h.update(chunk)
    return h.hexdigest()


def check_rows(rows):
    seen=set()
    for row in rows:
        name=row['path'];p=R/name
        assert name not in seen and not Path(name).is_absolute() and '..' not in Path(name).parts,name
        seen.add(name)
        assert p.is_file() and not p.is_symlink() and p.stat().st_size==row['bytes'] and digest(p)==row['sha256'],name
    return len(seen)


def verify():
    raw=json.loads((H/'development-raw-manifest.json').read_text())
    nraw=check_rows(raw['files'])
    freeze=json.loads((H/'development-freeze.json').read_text())
    for name,expected in freeze['code_sha256'].items():
        assert digest(H/'frozen-v1'/name)==expected,name
    state=json.loads((H/'development/state.json').read_text())
    result=json.loads((H/'adjudicated-summary.json').read_text())
    assert state['stage']=='completed' and result['complete'] and not result['missing']
    assert state['completed_pairs']==result['raw_pairs_replayed']==332
    assert len(result['rows'])==32 and not result['integer_oracle_native_mismatches']
    assert (result['adequate_workflows'],result['failed_workflows'],result['heldout_truth_supported'])==(1,9,10)
    assert not result['empirical_utility_comparison_supported'] and not result['strict_increment_over_all_strong_suites']
    assert result['logical_execution_reservation']==752
    assert result['logical_workflow_dispatches']==680
    assert result['actual_workflow_executions']==380 and result['exact_result_cache_hits']==300
    assert result['actual_calibration_executions']==284 and result['actual_decoder_executions']==380
    before=json.loads(subprocess.check_output(['git','show','94f63d0:'+str((H/'adjudicated-summary.json').relative_to(R))],cwd=R))
    assert result['confusion']==before['confusion'] and result['per_preparation']==before['per_preparation']
    assert result['per_policy_measured_union_runtime']==before['per_policy_measured_union_runtime']
    entry=json.loads((H/'portable-full-workflow-input.json').read_text())
    assert digest(R/entry['path'])==entry['sha256']
    native=json.loads((H/'fresh-full-workflow-qualification.json').read_text())
    assert native['input_sha256']==entry['sha256'] and native['passed'] and native['eligible_for_confirmation']
    assert native['native_nonblank_calls']==857 and native['calls'][0]['diagnostics']['blank_calls']==0
    assert native['native_blank_proxy_upper95']<=.1 and native['actual_workflow_executions']==2
    assert all(r['coverage']>=.99 for r in native['calls'])
    assert len(list((H/'full-workflow-control-calls').glob('call-*.json')))==4
    assert all(json.loads(p.read_text())['stage']=='completed' for p in (H/'full-workflow-control-calls').glob('call-*.json'))
    pins=json.loads((H/'full-workflow-code-pins.json').read_text())
    for name,expected in pins['code_sha256'].items():
        p=H/'frozen-v1'/name if name=='tool/src/codedcal/workflow_validation.py' else R/name
        assert digest(p)==expected,name
    decision=json.loads((H/'confirmation-decision.json').read_text())
    assert not decision['run_comparative_confirmation'] and decision['external_comparative_logical_calls']==0
    assert decision['conservative_total_charge']==84+2+1+752+4==843<=1600
    receipts=json.loads((H/'fullfield-context-receipts.json').read_text())
    assert receipts['requests']==32 and receipts['transport_bytes']==124246528
    context=json.loads((H/'native-context-verification.json').read_text())
    assert context['passed'] and len(context['records'])==8 and not context['whole_TIFF_hashes_recomputed_offline']
    # Imported science is verified separately by reproduce/verify.py. All other
    # files present at the preceding baseline are checked against their Git blobs.
    imported={r['path'] for r in json.loads((R/'provenance/import-manifest.json').read_text())['files']}
    navigation={r['path'] for r in json.loads((R/'provenance/repository-manifest.json').read_text())['files']}
    extra_count=0;extra_bytes=0
    tree=subprocess.check_output(['git','ls-tree','-r','-z',BASELINE],cwd=R)
    for record in tree.split(b'\0'):
        if not record:continue
        metadata,name=record.split(b'\t',1);mode,kind,oid=metadata.decode().split();name=name.decode()
        if name in imported or name in navigation:continue
        assert mode in ('100644','100755') and kind=='blob',name
        p=R/name;size=p.stat().st_size
        assert digest(p,'sha1',('blob '+str(size)+'\0').encode())==oid,name
        extra_count+=1;extra_bytes+=size
    assert (extra_count,extra_bytes)==(2357,446645341),(extra_count,extra_bytes)
    historical=json.loads((H/'completion-historical-closure.json').read_text())
    assert historical['passed'] and historical['imported_files']==19089
    assert historical['current_code_and_navigation_files']==len(navigation)==62
    original_cli=subprocess.check_output(['git','show',BASELINE+':tool/src/codedcal/cli.py'],cwd=R)
    assert original_cli==(R/'tool/src/codedcal/cli.py').read_bytes()
    manifest=json.loads((H/'completion-manifest.json').read_text())
    sealed=check_rows(manifest['files'])
    return dict(passed=True,raw_development_files_verified=nraw,sealed_files_verified=sealed,
        frozen_source_files_verified=len(freeze['code_sha256']),original_curated_files_verified=19089,
        preceding_code_and_navigation_files_verified=62,
        additional_preexisting_files_verified=extra_count,additional_preexisting_bytes_verified=extra_bytes,
        original_estimator_cli_byte_identical=True,analysis_classifications_and_runtimes_unchanged=True,
        complete_development=True,readout_aware_increment_supported=False,
        native_culture_workflow_qualified=True,independent_comparison_completed=False,
        external_comparative_logical_calls=0,conservative_total_charge=843,
        new_network_calls=0,new_fits_during_verification=0,physical_validation_added=False)


def seal():
    omitted={'completion-manifest.json','completion-verification.json'}
    paths=[p for p in H.rglob('*') if p.is_file() and p.name not in omitted and '__pycache__' not in p.parts and not p.name.endswith('.tmp')]
    paths += [R/'tool'/p for p in ['src/codedcal/workflow_validation.py','src/codedcal/workflow_cli.py',
        'src/codedcal/workflow_reporting.py','src/codedcal/workflow_bardensr.py','src/codedcal/workflow_decisions.py',
        'tests/test_workflow_validation.py','tests/test_workflow_reporting.py','docs/workflow_challenges.md']]
    rows=[dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(set(paths))]
    (H/'completion-manifest.json').write_text(json.dumps(dict(scope='Completed bounded study and additive software; immutable earlier evidence checked separately',files=rows),indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--seal',action='store_true');args=parser.parse_args()
    if args.seal:seal()
    report=verify()
    (H/'completion-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
