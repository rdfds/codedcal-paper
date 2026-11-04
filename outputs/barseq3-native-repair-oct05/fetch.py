"""One bounded, checksum-verified public file; never downloads other archives."""
import hashlib
import json
from pathlib import Path
import shutil
import requests
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
CACHE=ROOT.parent/'barseq3-inputs-oct05.nosync'
source=next(f for f in json.loads((HERE.parent/'independent-comparison-oct05/barseq3-zenodo.json').read_text())['body']['files'] if f['key'].startswith('adult_cb_raw'))
path=CACHE/'adult-cb-full.hdf5';expected='a03d84ad28824b9fef6ff2ffd514eee3'
assert source['size']==934285696 and source['checksum']=='md5:'+expected
assert shutil.disk_usage(CACHE).free>=3*1024**3
if path.exists():
    assert path.stat().st_size==source['size']
    h=hashlib.md5()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024**2),b''):h.update(b)
    assert h.hexdigest()==expected
    count=source['size']
else:
    temp=path.with_suffix('.partial');h=hashlib.md5();count=0
    with requests.get(source['links']['self'],params={'native_full_verified':'oct05'},stream=True,timeout=(20,60)) as r:
        r.raise_for_status();assert r.status_code==200,'partial response not accepted as whole source'
        with temp.open('wb') as f:
            for b in r.iter_content(1024**2):
                count+=len(b);assert count<=source['size'],'published size cap exceeded'
                h.update(b);f.write(b)
    assert count==source['size'] and h.hexdigest()==expected
    temp.replace(path)
(HERE/'full-input-verification.json').write_text(json.dumps(dict(source=source,bytes=count,md5=h.hexdigest(),whole_source_verified=True,path=str(path.resolve()),protocol_sha256=hashlib.sha256((HERE/'PROTOCOL.md').read_bytes()).hexdigest()),indent=2)+'\n')
print('Verified full source',count,h.hexdigest(),flush=True)
