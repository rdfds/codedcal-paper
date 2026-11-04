from pathlib import Path
import sys,json,time
import numpy as np
from scipy.ndimage import gaussian_filter
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path[:0]=[str(R/'tool/src'),str(R.parent/'bardensr-upstream.nosync')]
from codedcal.workflow_bardensr import make_adapter
from codedcal.workflow_validation import ImageCase,warp,canonical,pad_for_challenges,write_json
out=H/'native-numerical-controls.json';assert not out.exists();rows=[]
for colors in [2,4]:
    rng=np.random.default_rng(110126);n=8*colors;image=np.zeros((1,128,128),np.float32)
    for y,x in rng.integers(32,96,size=(12,2)):image[0,y,x]+=rng.uniform(1,3)
    image=gaussian_filter(image,(0,1,1));native=np.repeat(image[None],n,axis=0)
    imposed=np.zeros((n,3));imposed[:,1]=np.tile(np.where(np.arange(colors)<colors//2,3,-3),8)
    case=ImageCase(warp(native,imposed),np.ones_like(native),np.ones((n,1)),np.repeat(np.arange(8),colors),np.tile(np.arange(colors),8),np.array([0]),np.array([],int),'synthetic','synthetic')
    adapter=make_adapter();start=time.perf_counter();maps,diag=adapter.calibrate(pad_for_challenges(case),'group_search')
    rms=float(np.sqrt(np.mean(np.square(canonical(maps+imposed)).sum(axis=1))))
    rows.append(dict(colors=colors,residual_rms=rms,native_loss_checked=diag['native_loss_equivalence_checked'],passed=rms<1e-6,diagnostics=diag,wall_seconds=time.perf_counter()-start));print(colors,rms,flush=True)
write_json(out,dict(passed=all(r['passed'] for r in rows),records=rows,synthetic_only=True))
