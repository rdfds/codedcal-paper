from pathlib import Path
import json,sys,time,hashlib
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from native_adapters import make_adapter
from codedcal.workflow_validation import ImageCase,qualify_input,write_json
entries=json.loads((HERE/'exposed-inputs.json').read_text())['records']
start=time.monotonic();OUT=HERE/'native-qualification-v3';OUT.mkdir(exist_ok=True)
for e in entries:
    assert hashlib.sha256(Path(e['path']).read_bytes()).hexdigest()==e['sha256']
    case=ImageCase.load(e['path'])
    for decoder in (['correlation'] if e['name'].startswith('barseq3-') else ['correlation','regression']):
        path=OUT/(e['name']+'-'+decoder+'.json');assert not path.exists()
        row=qualify_input(case,make_adapter({'decoder':decoder}));row.update(input_sha256=e['sha256'],decoder=decoder)
        write_json(path,row);print(e['name'],decoder,row['passed'],row['native_nonblank_calls'],round(time.monotonic()-start,2),flush=True)
