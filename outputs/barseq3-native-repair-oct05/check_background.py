"""Check streaming port against actual dependency, in a short TF process."""
import os
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL','2')
os.environ.setdefault('TF_NUM_INTRAOP_THREADS','2')
os.environ.setdefault('TF_NUM_INTEROP_THREADS','1')
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
sys.path[:0]=[str(ROOT.parent/'bardensr-upstream.nosync'),str(HERE)]
import bardensr
from condition import CACHE,sha,write
checks=[]
with np.load(CACHE/'native-background-control.npz') as data:
    for dtype in ['32','64']:
        for sigma in [5,20]:
            native=bardensr.preprocessing.background_subtraction(data['a'+dtype],[0,sigma,sigma])
            port=data[f'b{dtype}s{sigma}'];error=float(np.max(np.abs(native-port)))
            assert error<(2e-6 if dtype=='32' else 1e-12)
            checks.append(dict(dtype=dtype,sigma=sigma,max_absolute_error=error))
write(HERE/'background-controls.json',dict(passed=True,checks=checks,
    upstream_background_sha256=sha(ROOT.parent/'bardensr-upstream.nosync/bardensr/preprocessing/preprocessing_tf.py'),
    condition_sha256=sha(HERE/'condition.py'),behavior='last-axis only, zero padding, native 3-sigma kernel'))
