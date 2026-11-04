"""Acquire only the immutable public source used by the audited AST adapter."""
import hashlib
import json
from pathlib import Path
import requests
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent
path=ROOT.parent/'barseq3-inputs-oct05.nosync/source-audit/barseq3-Snakefile_nopatch.txt'
url='https://raw.githubusercontent.com/kebschulllab/BARseq3/32df2c4793b8ab0534b73cc56a506ae16d4048a8/BARseq3_pipeline/gene_only_raw_image_pipeline/Snakefile_nopatch'
expected='31c8aa07c82e5e51968de0ca3d5111dd45c0849326b1c1c49b966b65ad9d5a94'
if path.exists():data=path.read_bytes()
else:
    response=requests.get(url,timeout=40);response.raise_for_status();data=response.content
    assert len(data)<512*1024
assert hashlib.sha256(data).hexdigest()==expected,'Pinned source mismatch; do not execute'
if not path.exists():
    path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
(HERE/'source-verification.json').write_text(json.dumps(dict(url=url,sha256=expected,bytes=len(data),
    commit='32df2c4793b8ab0534b73cc56a506ae16d4048a8',public_source_verified=True),indent=2)+'\n')
print('Verified pinned public native source',len(data),expected,flush=True)
