"""Finite-profile controls; no public pixels or comparative outcome labels."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.ndimage import shift
ROOT=Path(__file__).resolve().parent.parent.parent;sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal import output_risk as O
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run():
    book=np.zeros((8,16));
    for a in range(4):
        for b in range(4):book[a,4*a+b]=book[4+b,4*a+b]=1
    images=np.zeros((8,1,128,128));
    for y in range(20,108,8):
        for i,x in enumerate(range(20,108,4)):images[0,0,y,x]=1;images[4+i%4,0,y,x]=1
    support=np.ones_like(images,bool);imposed=np.zeros((8,3),int);imposed[4:,2]=4
    observed=np.stack([shift(a,t,order=0,mode='constant') for a,t in zip(images,imposed)])
    maps=np.stack([np.zeros((8,3),int),-imposed,imposed]);reg=((0,1),(40,88),(40,88));blocks=(1,12,12)
    profile=O.profile_losses(observed,support,book,maps,reg,blocks);comp=O.compatibility(profile);assert not comp['abstained'] and all(comp['compatible']);assert max(comp['losses'])<1e-6
    def catalog(data):
        z,y,x=np.where(data[0]>0);keep=(x>=40)&(x<88)&(y>=40)&(y<88);x,y,z=x[keep],y[keep],z[keep];values=data[:,z,y,x].reshape(2,4,-1);word=4*values[0].argmax(axis=0)+values[1].argmax(axis=0);return dict(xyz=np.c_[x,y,z].tolist(),words=word.tolist())
    cats=[catalog(np.stack([shift(a,t,order=0,mode='constant') for a,t in zip(observed,m)])) for m in maps];truth=catalog(images);assert O.catalog_distance(cats[0],truth)['identity_loss']==1 and O.catalog_distance(cats[1],truth)['identity_loss']==0
    risk=O.output_risk(cats,comp,0,[0]);assert risk['flagged'] and risk['identity_risk']==1 and risk['ordinary_multistart_identity_disagreement']==0
    # Pure global intensity gain leaves normalized residual and compatibility unchanged.
    scaled=O.compatibility(O.profile_losses(observed*7,support,book,maps,reg,blocks));assert np.allclose(comp['losses'],scaled['losses'],atol=1e-7) and comp['compatible']==scaled['compatible']
    shifted=maps+np.array([2,3,4]);assert np.array_equal(O.canonical_maps(maps),O.canonical_maps(shifted))
    assert O.compatibility(profile)==O.compatibility(profile)
    masked=O.profile_losses(observed,np.zeros_like(support),book,maps,reg,blocks);assert O.compatibility(masked)['abstained']
    unsupported=support.copy();unsupported[7,:,45:55,45:55]=False;partial=O.profile_losses(observed,unsupported,book,maps,reg,blocks);assert 0<partial['coverage']<1 and partial['residual'].shape[1]==16
    confident=O.compatibility(O.profile_losses(images,support,book,maps[:1],reg,blocks));native=O.output_risk([truth],confident,0,[0]);assert not native['flagged'] and not native['molecular_correctness_certified']
    corrupt=truth|dict(words=[(w+1)%4 for w in truth['words']]);invisible=O.catalog_distance(truth,corrupt);assert invisible['identity_loss']==1 and invisible['composition_tv']==0 and not native['flagged']
    try:O.canonical_maps(maps+.25)
    except ValueError:pass
    else:raise AssertionError('Fractional transforms accepted')
    return dict(passed=True,observable_false_consensus_identity_loss=1.,observable_false_consensus_flagged=True,ordinary_consensus_disagreement=0.,false_consensus_is_constructed_fit_bank_not_reexecution_of_all_methods=True,exact_inverse_identity_loss=0.,global_gain_invariant=True,first_frame_gauge_invariant=True,common_support_checked=True,block_resampling_deterministic=True,unsupported_abstention=True,invisible_label_permutation_identity_loss=invisible['identity_loss'],invisible_label_permutation_composition_tv=invisible['composition_tv'],invisible_label_permutation_flagged=native['flagged'],RNA_correctness_certification=False,script_sha256=sha(__file__),core_sha256=sha(ROOT/'tool/src/codedcal/output_risk.py'),protocol_sha256=sha(HERE/'PROTOCOL.md'))
if __name__=='__main__':
    out=run();(HERE/'controls.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
