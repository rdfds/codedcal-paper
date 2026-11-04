"""Six-round global-core and independent matching/oracle controls."""
import numpy as np
from scipy.ndimage import gaussian_filter
import adapter as A
import global_repair as G
import probe as P

def main():
    A.guard();codes,negative,*_=A.book();physical=codes[:,:998];rng=np.random.default_rng(61010)
    positions=np.array([(y,x) for y in range(29,100,10) for x in range(29,100,10)])[:60];ids=rng.choice(998,60,replace=False)
    images=np.zeros((24,128,128),np.float32)
    for (y,x),j in zip(positions,ids):images[:,y,x]+=physical[:,j]
    images=gaussian_filter(images,(0,1,1));results=[]
    for kind in ['zero','color']:
        imposed=P.shifts(kind);observed=A.R.warp(images,imposed);support=A.R.warp(np.ones_like(images),imposed)
        maps,diag=G.fit(observed,support,physical);exact=np.array_equal(A.R.warp(observed,maps)[:,16:112,16:112],images[:,16:112,16:112])
        r=dict(kind=kind,exact_input_recovery=exact,added_rms=float(np.sqrt(np.mean(np.sum((imposed+maps)**2,axis=1)))),diagnostics=diag);results.append(r)
        A.write(P.DEST/'controls.json',dict(passed=False,global_controls=results));print('SIX ROUND CONTROL',kind,exact,r['added_rms'],flush=True);assert exact
    a=dict(xyz=[[0,0,0],[1.5,0,0]],words=[0,1]);b=dict(xyz=[[.75,0,0],[-.9,0,0]],words=[1,0]);scores=P.metrics(a,b,1.)
    assert scores['identity_f1']==1 and scores['correct_labels']==2
    labels=dict(xyz=[[0,0,0],[1.5,0,0]],words=[1,0]);corruption=P.metrics(a,labels,.1)
    assert corruption['identity_f1']==0 and corruption['composition_tv']==0
    empty=P.metrics(dict(xyz=[],words=[]),dict(xyz=[],words=[]),1.);assert empty['identity_f1'] is None
    assert np.array_equal(np.sort(P.shifts('common'),axis=0),np.sort(P.shifts('color'),axis=0))
    A.write(P.DEST/'controls.json',dict(passed=True,global_controls=results,three_dimensional_maximum_cardinality_matching=True,label_corruption_identity_zero_composition_unchanged=True,empty_reference_undefined=True,complete_offset_distributions_equal=True,seed=61010,probe_sha256=A.sha(P.__file__),script_sha256=A.sha(__file__),development=True,independent_confirmation=False))
    print('PASSED six-round core, matching, corruption and distribution controls',flush=True)

if __name__=='__main__':main()
