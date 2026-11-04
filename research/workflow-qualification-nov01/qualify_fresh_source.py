from pathlib import Path
import sys,json,time,hashlib,signal
import numpy as np
from scipy.stats import chi2
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path[:0]=[str(R/'tool/src'),str(R.parent/'bardensr-upstream.nosync')]
from codedcal.workflow_validation import ImageCase,qualify_input,write_json
sys.path.insert(0,str(H))
from native_serval_adapter import make_adapter
entry=json.loads((H/'fresh-input-v2.json').read_text());p=Path(entry['path']);assert hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256'];case=ImageCase.load(p)
O=H/'fresh-qualification-source';O.mkdir(exist_ok=True)
signal.signal(signal.SIGALRM,lambda *args: (_ for _ in ()).throw(TimeoutError('Finite30-minute native-qualification limit')));signal.alarm(1800)
for decoder in ['correlation']:
    target=O/(decoder+'.json');assert not target.exists()
    row=qualify_input(case,make_adapter({'decoder':decoder}));native=row['calls'][0];diag=native.get('diagnostics',{});b=diag.get('blank_calls');n=diag.get('nonblank_calls');nblank=len(case.blank_classes);nclasses=len(np.unique(case.class_ids))
    upper=float(.5*chi2.ppf(.95,2*(b+1))*(nclasses-nblank)/nblank/n) if nblank and n and b is not None else None
    row.update(decoder=decoder,input_sha256=entry['sha256'],native_blank_proxy_upper95=upper,blank_proxy_assumptions='Named blank-codewords as a proxy; not molecular ground truth',eligible_for_confirmation=row['passed'] and (upper is not None and upper<=.1 if nblank else True))
    write_json(target,row);print(decoder,row['passed'],row['native_nonblank_calls'],upper,row['eligible_for_confirmation'],flush=True)
