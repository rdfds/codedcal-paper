"""Known-gain/source-label control of the actual CoSiD implementation."""
from pathlib import Path
import sys,hashlib,json,random
import numpy as np
import pandas as pd
H=Path(__file__).resolve().parent;R=H.parents[1];sys.path[:0]=[str(R/'tool/src'),str(R.parent/'serval-upstream.nosync')]
from codedcal.workflow_validation import write_json
from serval.codebook import Codebook
from serval.image import ImageStack
from serval.decode.pixel import NearestNeigbourPixelDecoder,CosineOptimizedPixelDecoder
out=H/'native-source-gain-control.json';assert not out.exists()
rng=np.random.default_rng(110126);words=[]
while len(words)<32:
 w=tuple(rng.integers(0,4,4))
 if w not in words:words.append(w)
b=np.zeros((16,32));a=np.zeros((16,128,128));truth=[]
for j,w in enumerate(words):b[np.arange(4)*4+np.array(w),j]=1
true_gain=rng.uniform(.85,1.15,16)
for i in range(96):
 j=i%32;y=10+10*(i//10);x=10+10*(i%10);a[:,y:y+2,x:x+4]=20*(b[:,j]*true_gain)[:,None,None];truth.append((j,y+.5,x+1.5))
book=Codebook(pd.DataFrame(b.T));decoder=CosineOptimizedPixelDecoder(book,NearestNeigbourPixelDecoder(book,max_dist=.5167,min_norm=0),init_scaling_factors=np.ones(16),min_area=5,max_area=12,max_size=100000,penalty_entropy=0,penalty_l2=0)
random.seed(110126)
for _ in range(10):decoder.fit([ImageStack(a)])
ratio=decoder.scaling_factors/decoder.scaling_factors.mean();target=true_gain/true_gain.mean();error=float(np.max(np.abs(ratio-target)))
decoder.decoder=NearestNeigbourPixelDecoder(book,max_dist=.5167,min_norm=1);table=decoder.predict(ImageStack(a)).spots
observed=sorted(map(tuple,table[['barcode_id','y','x']].to_numpy(float)));exact=observed==sorted(truth)
write_json(out,dict(passed=error<=.03 and exact,max_relative_gain_error=error,known_source_labels_and_centroids_exact=exact,source_objects=len(truth),recovered_objects=len(table),true_relative_gains=target,estimated_relative_gains=ratio,seed=110126,synthetic_only=True,decoder_only_executions=1,source_code_sha256=hashlib.sha256((R.parent/'serval-upstream.nosync/serval/decode/pixel/cosine.py').read_bytes()).hexdigest()))
print('NATIVE GAIN CONTROL',error,exact,len(table),flush=True)
