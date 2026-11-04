"""Known operator/endpoint controls; no public images, no TensorFlow inference."""
import ast
from pathlib import Path
import sys
import numpy as np
from scipy.ndimage import shift,gaussian_filter
from scipy.signal import fftconvolve
from condition import HERE,ROOT,write,sha
sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal.evaluation import spatial_match
ns=dict(np=np,shift=shift,fftconvolve=fftconvolve,spatial_match=spatial_match)
tree=ast.parse((HERE/'recovery.py').read_text())
for name in ['offsets','warp','bounded','refine','metrics']:
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
    exec(compile(ast.Module(body=[node],type_ignores=[]),'frozen recovery control','exec'),ns)
rng=np.random.default_rng(505);image=np.zeros((128,128))
for x,y in rng.integers(32,96,size=(20,2)):image[y,x]+=rng.uniform(.5,1.5)
image=gaussian_filter(image,1.4);imposed=np.array([4,2]);observed=ns['warp'](image[None],imposed[None])[0]
correction=ns['bounded'](image,observed)
assert np.max(np.abs(correction+imposed))<.1
pixel_checks={}
for kind in ['common','color']:
    offsets=ns['offsets'](kind);raw=np.repeat(image[None],40,axis=0)
    recovered=ns['warp'](ns['warp'](raw,offsets),-offsets)
    assert np.array_equal(recovered[:,16:112,16:112],raw[:,16:112,16:112])
    pixel_checks[kind]=dict(exact_oracle_decoder_window=True,covariance=(offsets.T@offsets/40).tolist())
xy=np.array([[40,40],[80,80]],float);labels=np.array([0,1]);names=['A','B'];ids=np.array([0,1])
identity=ns['metrics'](xy,labels,xy,labels,ids,names,5)
swap=ns['metrics'](xy,labels,xy,labels[::-1],ids,names,5)
assert identity['gene_class_identity_f1']==1 and swap['gene_class_identity_f1']==0
assert swap['spatial_f1']==1 and swap['composition_tv']==0
synonym=ns['metrics'](xy[:1],np.array([0]),xy[:1],np.array([1]),np.array([0,0,1]),['A','B'],5)
assert synonym['gene_class_identity_f1']==1 and synonym['sequence_identity_f1']==0
write(HERE/'recovery-operator-controls.json',dict(passed=True,recovery_sha256=sha(HERE/'recovery.py'),
    synthetic_source_only=True,known_integer_translation_max_error=float(np.max(np.abs(correction+imposed))),
    pixel_controls=pixel_checks,label_corruption_detected=True,spatial_and_composition_invariant=True,synonymous_barcodes_grouped=True))
print('Recovery operator and endpoint controls passed',flush=True)
