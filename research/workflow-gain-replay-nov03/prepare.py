"""Freeze all existing inputs, records, inference code and decisions before replay."""
from pathlib import Path
import datetime
import importlib.metadata
import json
import platform
import subprocess
import sys
from replay_common import *


def main():
    assert not (H/'freeze.json').exists()
    original=load(P/'freeze.json')
    assert load(P/'run/state.json')['stage']=='completed'
    seals={}
    for name in ['workflow-validation-completion-nov03/completion-manifest.json',
                 'workflow-failure-audit-nov03/posthoc-manifest.json',
                 'workflow-failure-audit-nov03/native-completion-manifest.json']:
        path=R/'research'/name
        seals[str(path.relative_to(R))]=dict(sha256=sha(path),files_verified=verify_manifest(path))
    for name,digest in original['code_sha256'].items():assert sha(R/name)==digest,name
    assert sha(R/original['input']['path'])==original['input']['sha256']
    assert load(H/'software-control-receipt.json')['actual_native_decoder_predictions']==4
    calls=[(p,load(p)) for p in sorted((P/'run/calls').glob('*.json'))]
    case=pad_for_challenges(ImageCase.load(R/original['input']['path']))
    entries=[]
    for policy in ['barcode_model','group_search']:
        native=load(P/'run/trials'/policy/'native.json')
        gate=load(P/'run/qualification'/(policy+'.json'))
        assert gate['primary_native_eligible']
        challenges=[None]+[c for c in original['challenges'] if c['id'] not in original['fractional']]
        assert len(challenges)==13
        for ch in challenges:
            name='native' if ch is None else ch['id']
            path=P/'run/trials'/policy/(name+'.json');row=load(path)
            old=row if ch is None else row['runs'][0]
            if ch is not None:
                assert row['integer_oracle_native_exact'] and row['runs'][1]['catalog']==native['catalog']
            raw=observed(case,ch);key=cache_key(raw,original['adapter_config'],policy)
            matches=[(p,c) for p,c in calls if c['stage']=='completed' and c['policy']==policy and c['input_fingerprint']==key]
            assert matches,(policy,name)
            for p,c in matches:
                assert c['catalog']==old['catalog'] and c['corrections']==old['corrections'] and c['coverage']==old['coverage']
                assert c['diagnostics']['fitted_source_gains']==old['diagnostics']['fitted_source_gains']
            entries.append(dict(id=policy+'/'+name,policy=policy,challenge=ch,input_fingerprint=key,
                source_path=str(path.relative_to(R)),source_sha256=sha(path),
                matching_completed_calls=[dict(path=str(p.relative_to(R)),sha256=sha(p)) for p,c in matches]))
            del raw
    env={k:importlib.metadata.version(k) for k in list(original['environment'])+['scikit-learn']}
    assert all(env[k]==v for k,v in original['environment'].items())
    own=sorted(H.glob('*.py'))+[H/'PROTOCOL.md',H/'software-control-receipt.json',H/'software-tests.log',H/'software-tests-repaired.log']
    codes={str(p.relative_to(R)):sha(p) for p in own}
    for name in ['tool/src/codedcal/workflow_validation.py','tool/src/codedcal/workflow_gain.py','tool/src/codedcal/workflow_decisions.py',
                 'research/workflow-validation-controller-nov02/common.py','research/workflow-memory-repair-nov02/supervise_finite.py']:
        codes[name]=sha(R/name)
    # Include the entire existing upstream Python source pins, not just the wrapper.
    for name,digest in original['serval_sha256'].items():assert sha(R.parent/'serval-upstream.nosync'/name)==digest,name
    freeze=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        user_authorization='ok do it and do it well: one fixed-fit replay of all24 integer cases plus two native controls',
        preceding_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),
        original_freeze_path=str((P/'freeze.json').relative_to(R)),original_freeze_sha256=sha(P/'freeze.json'),
        input=original['input'],adapter_config=original['adapter_config'],environment=env,python=platform.python_version(),
        code_sha256=codes,serval_sha256=original['serval_sha256'],preceding_seals=seals,
        cases=entries,policies=['barcode_model','group_search'],strategies=original['strategies'],heldout=original['heldout'],
        fractional=original['fractional'],fractional_exclusion=original['fractional_exclusion'],
        originally_assigned_policies=original['policies'],preparation=original['preparation'],
        original_reference_catalogs=True,coordinate_engineering_tolerance=1e-12,
        criteria=dict(native_minimum=50,reference_minimum=10,coverage_minimum=.99,blank_proxy_upper95_maximum=.1,
            identity_and_composition_tolerance=.05,matching_radius=3.,sensitivity_radii=[1.5,5.]),
        maximum_decoder_predictions=52,software_decoder_predictions=4,prior_conservative_reservation=1159,
        total_conservative_reservation=1215,maximum_total=1600,maximum_wall_seconds=1800,
        minimum_start_free_bytes=2*1024**3,minimum_running_free_bytes=512*1024**2,maximum_retained_bytes=64*1024**2,
        utility_minimum_adequate=5,utility_minimum_failed=5,utility_minimum_preparations=2,
        comparative_outcomes_exposed=True,independent_confirmation=False,physical_validation=False)
    assert len(entries)==26 and 1159+4+2*len(entries)==1215<=1600
    write_json(H/'freeze.json',freeze)
    print(json.dumps(dict(cases=len(entries),planned_decoder_predictions=52,environment=env,preceding_seals=seals),indent=2))


if __name__=='__main__':main()
