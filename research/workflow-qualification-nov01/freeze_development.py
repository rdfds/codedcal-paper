from pathlib import Path
import sys,json,hashlib,subprocess
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path.insert(0,str(R/'tool/src'))
from codedcal.workflow_validation import ImageCase,challenge_plan,write_json
out=H/'development-freeze.json';assert not out.exists()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
entries=json.loads((H/'exposed-inputs.json').read_text())['records'];cases=[]
for e in entries:
    assert sha(e['path'])==e['sha256'];case=ImageCase.load(e['path']);decoders=['correlation'] if e['name'].startswith('barseq3-') else ['correlation','regression']
    gates={d:json.loads((H/'native-qualification-v3'/(e['name']+'-'+d+'.json')).read_text()) for d in decoders}
    cases.append(dict(name=e['name'],preparation=e['preparation'],path=e['path'],sha256=e['sha256'],decoders=decoders,
        qualified={d:gates[d]['passed'] for d in decoders},qualification_sha256={d:sha(H/'native-qualification-v3'/(e['name']+'-'+d+'.json')) for d in decoders},plan=challenge_plan(case)))
paths=[R/'tool/src/codedcal'/p for p in ['workflow_validation.py','workflow_bardensr.py','workflow_decisions.py']]
paths += [H/p for p in ['native_adapters.py','run_development.py','summarize.py','PROTOCOL.md','QUALIFICATION_REVISION.md','SUPPORT_MEASURES.md']]
paths += list((R/'outputs/starfinder-public-confirmation-oct06').glob('*.py'))+list((R/'outputs/barseq3-native-repair-oct05').glob('*.py'))
freeze=dict(question='Qualification of computational catalog preservation at four challenge pairs per strategy',cases=cases,
 policies=['round','frame_template','barcode_model','group_search'],seed=110126,identity_tolerance=.05,composition_tolerance=.05,matching_radius=3.,minimum_native=50,minimum_reference=10,minimum_coverage=.99,
 logical_limit=1600,planned_development_logical=1128,worker_seconds_limit=10800,decoder_regression_iterations=100,registration_iterations=50,
 calibration_starts=3,thread_limit=4,known_native_regression_compute_exclusions=['barseq3-0/regression','barseq3-1/regression'],
 code_sha256={str(p.relative_to(R)):sha(p) for p in sorted(set(paths))},git_parent=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),
 no_comparative_outcomes_examined=True,statistical_scope='Descriptive within-preparation workflow decisions; fields/decoders/challenges are not independent biological replicates',
 budget_scope='Equal maximum full-workflow invocations; measured actual resource/time and exact cache reuse reported separately, not claimed identical wall-clock seconds',
 threshold_learning='LOPO scalar QC/sensitivity thresholds: train empirical false alerts<=.1, require>=5 adequate training workflows, otherwise abstain; no held-out-preparation training')
freeze['native_code_sha256']={str(p.relative_to(R.parent/'bardensr-upstream.nosync')):sha(p) for p in sorted((R.parent/'bardensr-upstream.nosync/bardensr').rglob('*.py'))}
freeze['qualified_workflows']=sum(sum(c['qualified'].values()) for c in cases)*4
assert json.loads((H/'native-numerical-controls.json').read_text())['passed']
write_json(out,freeze);print('FROZEN',len(cases),'inputs',sum(len(c['decoders']) for c in cases)*4,'workflows',freeze['planned_development_logical'])
