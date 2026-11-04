"""Qualification before added conventional-comparator native outcomes."""
import numpy as np
from scipy.ndimage import gaussian_filter
import adapter as A
import local_trace as L
import probe as P
import source_comparators as S

def main():
    L.guard();rng=np.random.default_rng(61011);a=np.zeros((30,128,128),np.float32)
    for z,y,x in [(6,35,35),(14,63,63),(22,91,91)]:a[z,y,x]=rng.uniform(.5,1)
    a=gaussian_filter(a,(.7,1,1));images=np.tile(a[None],(24,1,1,1));support=np.ones_like(images,np.uint8);rows=[]
    for kind in ['zero','common','color']:
        imposed=P.shifts(kind);obs=P.warp(images,imposed)
        for roundwise in [True,False]:
            if roundwise and kind=='color':continue
            maps=S.global_maps(obs,roundwise);out=S.move(obs,maps);exact=np.array_equal(out[:,:,16:112,16:112],images[:,:,16:112,16:112]);assert exact
            rows.append(dict(kind=kind,roundwise=roundwise,exact_interior_recovery=True,maps_zyx=maps.tolist()))
    for policy in S.POLICIES:
        out,alpha,diag=S.apply(images,support,policy);assert np.array_equal(out,images) and np.array_equal(alpha,support)
    A.write(S.DEST/'controls.json',dict(passed=True,geometry_controls=rows,zero_local_identity_all_policies=True,prior_qualified_pull_sign_and_demons_translation=True,script_sha256=A.sha(S.__file__),control_script_sha256=A.sha(__file__),protocol_sha256=A.sha(A.HERE/'SOURCE_COMPARATOR_PROTOCOL.md'),source_files={str(p.relative_to(A.ROOT)):A.sha(p) for p in [A.AUDIT/'starfinder-author-source/code-base/src/RegisterImagesGlobal.m',A.AUDIT/'starfinder-author-source/code-base/src/RegisterImagesLocal.m']},source_not_exact_matlab=True,seed=61011))
    print('PASSED source-round and 3D frame-template geometry controls',flush=True)
if __name__=='__main__':main()
