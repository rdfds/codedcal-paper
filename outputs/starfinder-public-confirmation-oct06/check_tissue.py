"""Four-round complete-book eligibility controls before tissue pixel exposure."""
import json
import numpy as np
from scipy.ndimage import gaussian_filter
import adapter as A
import tissue as T
import probe as P
import global_repair as G

def main():
    codes,words=T.book();rng=np.random.default_rng(61012);genes=rng.choice(64,60,replace=False)
    centers=np.array([(y,x) for y in range(29,100,10) for x in range(29,100,10)])[:60];depth=np.array([6,14,22])[np.arange(60)%3]
    images=np.zeros((16,30,128,128),np.float32)
    for j,(y,x),z in zip(genes,centers,depth):images[:,z,y,x]+=codes[:,j]
    images=gaussian_filter(images,(0,.7,1,1));images/=images.max();support=np.ones_like(images,np.uint8)
    lookup={(int(x),int(y),int(z)):int(j) for j,(y,x),z in zip(genes,centers,depth)};results=[]
    for dec in A.DECODERS:
        result=T.catalog(images,support,codes,dec);correct=sum(lookup.get(tuple(p))==j for p,j in zip(result['xyz'],result['words']));f1=2*correct/(60+result['coding']);assert f1>=.99 and result['coding']==60 and result['unsupported_candidates']==0;results.append(dict(decoder=dec,f1=f1,coding=result['coding'],coverage=result['coverage'],runtime=result['runtime']))
    ambiguous=np.zeros_like(images)
    for (y,x),z in zip(centers,depth):ambiguous[:,z,y,x]=1
    ambiguous=gaussian_filter(ambiguous,(0,.7,1,1));ambiguous/=ambiguous.max()
    for dec in A.DECODERS:assert T.catalog(ambiguous,support,codes,dec)['coding']==0
    # A count-preserving identity permutation has unchanged composition.
    reference=dict(xyz=[[0,0,0],[10,0,0]],words=[0,1]);corrupt=dict(xyz=reference['xyz'],words=[1,0]);metric=P.metrics(reference,corrupt,1.);assert metric['identity_f1']==0 and metric['composition_tv']==0
    plane=images[:,14];geometry=[]
    for kind in ['zero','common','color']:
        imposed=T.offsets(kind);observed=A.R.warp(plane,imposed);maps,diagnostics=G.fit(observed,A.R.warp(np.ones_like(plane),imposed),codes);exact=np.array_equal(A.R.warp(observed,maps)[:,16:112,16:112],plane[:,16:112,16:112]);assert exact;geometry.append(dict(kind=kind,exact_input_recovery=True,diagnostics=diagnostics))
    a=np.zeros((30,128,128),np.float32);a[6,35,35]=1;a[14,64,64]=.8;a[22,91,91]=.5;a=gaussian_filter(a,(.7,1,1));shared=np.tile(a[None],(16,1,1,1))
    for kind in ['zero','common','color']:
        imposed=T.offsets(kind);observed=P.warp(shared,imposed);maps=T.maps3d(observed);assert np.array_equal(T.move3d(observed,maps)[:,:,16:112,16:112],shared[:,:,16:112,16:112])
    out,alpha,diag=T.source_prepare(shared,np.ones_like(shared,np.uint8));assert np.array_equal(out,shared) and np.array_equal(alpha,np.ones_like(shared))
    assert sorted(map(tuple,T.offsets('common')))==sorted(map(tuple,T.offsets('color')))
    assert not np.any(np.bitwise_xor.reduce(words,axis=1)) and codes.shape==(16,64)
    result=dict(passed=True,synthetic=results,global_core_geometry=geometry,threeD_template_controls_passed=True,source_identity_warp_passed=True,all_color_ambiguity_abstained=True,count_preserving_corruption_detected=True,complete_distribution_matched=True,no_comparable_absent_codes=True,tissue_pixels_unopened=True,seed=61012,script_sha256=A.sha(T.__file__),control_script_sha256=A.sha(__file__))
    A.write(T.DEST/'controls.json',result)
    paths=[A.HERE/n for n in ['tissue.py','check_tissue.py','TISSUE_STUDY_PROTOCOL.md','initialization.py','probe.py','local_trace.py','volume_native.py','adapter.py']]+[T.DEST/'controls.json',A.HERE/'freeze.json',A.AUDIT/'starfinder-color-codebooks.json',A.HERE.parent/'barseq3-native-repair-oct05/global_repair.py']
    A.write(T.DEST/'freeze.json',dict(files={str(p.relative_to(A.ROOT)):A.sha(p) for p in paths},tissue_pixels_unopened=True,scope='conditional native catalog; no RNA truth or absent-code proxy',full_comparative_workflow_frozen=True,protocol_sha256=A.sha(A.HERE/'TISSUE_STUDY_PROTOCOL.md')))
    print('PASSED four-round complete64book, actual kernels, geometry, ambiguity controls',flush=True)
if __name__=='__main__':main()
