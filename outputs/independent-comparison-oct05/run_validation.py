"""Guard validation with complete development and a committed preflight freeze."""
import run as study
import subprocess

freeze=study.HERE/'preflight-freeze.json'
state=study.json.loads(freeze.read_text())
assert state['validation_comparison_qualified'], 'oracle/interpolation gate failed'
assert state['settings_sha256']==study.digest(study.HERE/'settings-freeze.json')
assert state['strong_amendment_sha256']==study.digest(study.HERE/'STRONG_COMMON_AMENDMENT.md')
relative=str(freeze.relative_to(study.ROOT))
tracked=subprocess.run(['git','show','HEAD:'+relative],cwd=study.ROOT,
                        check=True,capture_output=True).stdout
assert tracked==freeze.read_bytes(), 'preflight must be committed before validation'
study.run('validation')
