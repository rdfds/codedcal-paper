"""Independent NumPy evaluation of actual native correlation references."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','2')
import numpy as np
from scipy.ndimage import maximum_filter
from condition import HERE,CACHE,read,write,sha

classes=read(HERE/'codebook-classes.json')
B=np.zeros((40,len(classes['sequences'])),np.float32)
for j,s in enumerate(classes['sequences']):
    for r,c in enumerate(s):B[4*r+'GTCA'.index(c),j]=1
B/=np.sqrt(np.sum(B**2,axis=0,keepdims=True))
records=[]
for index in [0,1]:
    ref=read(HERE/f'recovery/field-{index}/reference-correlation.json')
    arrays=[]
    for r in range(10):
        p=CACHE/f'native-repair-source-{index}-round-{r}.npz';marker=read(HERE/f'field-{index}/round-{r}-demons.json')
        assert sha(p)==marker['output_sha256']
        with np.load(p) as x:arrays.extend(x['images'])
    X=np.stack(arrays)[:,16:112,16:112].astype(np.float32).reshape(40,-1)
    X/=np.sqrt(np.sum(X**2,axis=0,keepdims=True))+np.float32(.05)
    density=(X.T@B).reshape(96,96,-1)
    maxima=maximum_filter(density,size=(3,3,1),mode='constant',cval=np.finfo(np.float32).min)
    y,x,j=np.where((density==maxima)&(maxima>np.float32(.775)));x+=16;y+=16
    interior=(x>=24)&(x<104)&(y>=24)&(y<104);x,y,j=x[interior],y[interior],j[interior]
    order=np.lexsort((x,y,j));xy=np.c_[x[order],y[order]].astype(float).tolist();words=j[order].tolist()
    equal=xy==ref['xy'] and words==ref['words']
    records.append(dict(field=index,independent_numpy_calls=len(words),native_tf_calls=len(ref['words']),exact_calls_match=equal,
        native_reference_sha256=sha(HERE/f'recovery/field-{index}/reference-correlation.json')))
    assert equal,'Independent correlation calls differ; inspect numerical boundary before any claim'
write(HERE/'independent-correlation-verification.json',dict(fields=records,all_pass=True,script_sha256=sha(__file__),
    scope='Independent numerical kernel and peak recount; not independent biological reference'))
print(records,flush=True)
