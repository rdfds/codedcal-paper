"""Qualify observable-only sampling; reuse only identical native kernel tensors."""
import ast
import json
import numpy as np
from scipy.ndimage import gaussian_filter
import active_plane as P
import adapter as A

def main():
    A.guard();codes,negative,names,words=A.book();rng=np.random.default_rng(61006)
    rng.integers(4,196,(3,13,17),dtype=np.uint8)
    genes=rng.choice(998,60,replace=False)
    centers=np.array([(y,x) for y in range(29,100,10) for x in range(29,100,10)])[:60]
    def image(z):
        images=np.zeros((24,30,128,128),np.float32)
        for gene,(y,x) in zip(genes,centers):
            for r,c in enumerate(words[gene]):images[r*4+c,z,y,x]=1
        images=gaussian_filter(images,sigma=(0,.7,1.,1.));images/=images.max();return images
    original=image(14);moved=image(6);support=np.ones_like(moved,np.uint8)
    selected,diag=P.select(moved,support);assert selected==6
    restored=P.center(moved,selected)
    assert np.array_equal(restored[:,13:16,16:112,16:112],original[:,13:16,16:112,16:112]),'Kernel input differs'
    source_xy,source_calls,_=A.source_decode(restored,codes)
    ref_xy,ref_calls,_=A.source_decode(original,codes)
    assert np.array_equal(source_xy,ref_xy) and np.array_equal(source_calls,ref_calls) and len(source_calls)==60
    # Selection must ignore later-round corruption and all code/identity labels.
    altered=moved.copy();altered[4:]=np.flip(altered[4:],axis=1)
    assert P.select(altered,support)[0]==selected
    # Empty windows remain explicit; no favorable replacement window is used.
    z,empty=P.select(np.zeros_like(moved),support);assert z==4 and empty['selected_reference_spots']==0
    supported=support.copy();supported[:,5:8]=0
    z,masked=P.select(moved,supported);assert z!=6 and 6 not in masked['safe_planes']
    prior=json.loads((A.HERE/'controls.json').read_text());assert prior['passed']
    reuse=[]
    for decoder in ['correlation','regression']:
        record=next(r for r in prior['synthetic'] if r['decoder']==decoder)
        assert record['f1']==1 and record['negative']==0
        reuse.append(record|dict(reuse_from=str((A.HERE/'controls.json').relative_to(A.ROOT)),exact_decoder_input_equality=True,actual_execution_kernel_s=0))
    result=dict(passed=True,development=True,selected_z=selected,reference_spots=diag['selected_reference_spots'],central_tensor_exact=True,source_calls_exact=60,later_round_invariance=True,empty_and_support_controls=True,kernel_control_reuse=reuse,sampling_sha256=A.sha(P.__file__),original_adapter_sha256=A.sha(A.HERE/'adapter.py'),script_sha256=A.sha(__file__),independent_confirmation=False)
    A.write(P.DEST/'controls.json',result)
    paths=[A.HERE/'active_plane.py',A.HERE/'check_active_plane.py',A.HERE/'ACTIVE_PLANE_PROTOCOL.md',P.DEST/'controls.json',A.HERE/'freeze.json',A.HERE/'preparation.json']
    files=json.loads((A.HERE/'freeze.json').read_text())['files'].copy()
    files.update({str(path.relative_to(A.ROOT)):A.sha(path) for path in paths})
    A.write(P.DEST/'freeze.json',dict(files=files,core_and_decoder_unchanged=True,development_pixels_already_exposed=True,independent_test_frozen=False,source_field='Position351',selection_uses_first_round_only=True))
    print('PASSED active sampling controls; exact decoder-input reuse; not independent confirmation',flush=True)

if __name__=='__main__':main()
