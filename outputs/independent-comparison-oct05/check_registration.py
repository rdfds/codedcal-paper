"""Known-source convention control; never reads real benchmark arrays."""
import run as study
import numpy as np
from scipy.ndimage import gaussian_filter

rng=np.random.default_rng(190)
book=np.zeros((7,4,32),np.float64)
for r in range(7):
    book[r,rng.integers(0,4,32),np.arange(32)]=1
sources=np.zeros((32,200,200),np.float64)
for j in range(32):
    for _ in range(4):
        y,x=rng.integers(30,170,2)
        sources[j,y,x]=rng.uniform(1,3)
sources=gaussian_filter(sources,(0,1.2,1.2))
book=book.reshape(28,32)
x=np.einsum('fj,jyx->fyx',book,sources)+.00001
imposed=study.offsets('color',2,501)
observed=study.warp(x,imposed)
common=study.phase(observed)
template,steps=study.template_alignment(observed,common)
model,losses=study.code_model(observed,book,template)
error=lambda t:float(np.sqrt(np.mean(np.sum((t+imposed)**2,axis=1))))
result=dict(kind='simulation_truth_convention_control',
            imposed_rms=error(np.zeros_like(imposed)),
            common_rms=error(common),template_rms=error(template),
            model_rms=error(model),wrong_sign_model_rms=error(-model),
            oracle_rms=error(-imposed),model_losses=losses,
            applied_model_yx=model.tolist(),imposed_yx=imposed.tolist())
study.write(study.HERE/'registration-control.json',result)
assert result['model_rms'] < result['wrong_sign_model_rms'], 'model sign not qualified'
assert result['template_rms'] < .25, 'known-source template recovery failed'
assert result['model_rms'] < .5, 'known-source model recovery failed'
print({k:v for k,v in result.items() if k.endswith('rms')},flush=True)
