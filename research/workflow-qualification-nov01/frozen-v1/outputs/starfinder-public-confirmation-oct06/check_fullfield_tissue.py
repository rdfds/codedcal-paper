"""Whole-field coarse/fine sign and phase-grid controls before repair outcomes."""
import numpy as np
from scipy.ndimage import gaussian_filter
import adapter as A
import fullfield_tissue as F

def main():
    rng=np.random.default_rng(61013);a=np.zeros((30,128,128),np.float32)
    for _ in range(30):z,y,x=rng.integers([5,20,20],[25,108,108]);a[z,y,x]+=rng.uniform(.5,1)
    a=gaussian_filter(a,(.7,1,1));a=np.rint(a/a.max()*255).astype(np.uint8);rows=[]
    for imposed in [(0,0,0),(0,4,2),(0,3,1),(1,-5,3),(-2,7,-3)]:
        moving=A.integer_warp(a,imposed);maps,diag=F.source_map(a,moving);exact=np.array_equal(A.integer_warp(moving,maps)[4:26,16:112,16:112],a[4:26,16:112,16:112]);assert exact;rows.append(dict(imposed_zyx=list(imposed),selected_zyx=maps.tolist(),exact=True,diagnostics=diag))
    A.write(F.DEST/'controls.json',dict(passed=True,records=rows,script_sha256=A.sha(F.__file__),control_script_sha256=A.sha(__file__),protocol_sha256=A.sha(A.HERE/'FULLFIELD_TISSUE_PROTOCOL.md'),unchanged_complete_book_controls_sha256=A.sha(F.T.DEST/'controls.json'),tissue_already_exposed_development=True))
    print('PASSED whole-field even/odd/Z coarse-fine controls',flush=True)
if __name__=='__main__':main()
