"""Independent record/count/assignment and finite-profile bootstrap recount."""
import collections,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent

def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def comparison(a,b):
    ax=np.array(a['xyz']).reshape(-1,3);bx=np.array(b['xyz']).reshape(-1,3);aw=np.array(a['words']);bw=np.array(b['words']);assert len(ax)==len(aw) and len(bx)==len(bw)
    if not len(aw) and not len(bw):return 0.,0.
    if not len(aw) or not len(bw):return 1.,1.
    d=cdist(ax,bx);valid=d<=3;maximum=int((maximum_bipartite_matching(csr_matrix(valid.astype(int)),perm_type='column')>=0).sum());reward=(min(len(ax),len(bx))+1)*3;cost=np.zeros((len(ax)+len(bx),len(ax)+len(bx)));cost[:len(ax),:len(bx)]=np.where(valid,d-reward,reward);i,j=linear_sum_assignment(cost);ok=(i<len(ax))&(j<len(bx));i,j=i[ok],j[ok];ok=valid[i,j];i,j=i[ok],j[ok];assert len(i)==maximum;loss=1-2*int((aw[i]==bw[j]).sum())/(len(aw)+len(bw));negative=set(a.get('negative_ids',[]))|set(b.get('negative_ids',[]));ca=collections.Counter(int(w) for w in aw if w not in negative);cb=collections.Counter(int(w) for w in bw if w not in negative);na=sum(ca.values());nb=sum(cb.values());tv=sum(abs(ca[w]/na-cb[w]/nb) for w in ca.keys()|cb.keys())/2 if na and nb else (0. if not na and not nb else 1.);return loss,tv

def main():
    freeze=read(HERE/'freeze.json')
    for rel,digest in freeze['files'].items():assert sha(ROOT/rel)==digest,rel
    evals=list((HERE/'cases').glob('*/evaluation.json'));assert len(evals)==13
    verified=[];catalogs_checked=0;profiles_checked=0;prediction_checks=0
    for p in sorted(evals):
        prof=read(p.parent/'profile.json');c=prof['compatibility'];res=np.array(prof['residual']);energy=np.array(prof['energy']);counts=np.array(prof['block_counts']);assert res.shape==energy.shape and len(res)==len(prof['maps_zyx']);assert np.isfinite(res).all() and (res>=0).all() and (energy>=0).all()
        if not c['abstained']:
            use=(counts>0)&(energy.max(axis=0)>0);res=res[:,use];energy=energy[:,use];score=res.sum(axis=1)/np.maximum(energy.sum(axis=1),1e-20);assert np.allclose(score,c['losses'],atol=1e-12);best=int(np.argmin(score));rng=np.random.default_rng(61021);indices=rng.integers(0,int(use.sum()),(500,int(use.sum())));boot=np.stack([res[:,i].sum(axis=1)/np.maximum(energy[:,i].sum(axis=1),1e-20) for i in indices],axis=1);lower=np.quantile(boot-boot[best],.025,axis=1);assert np.array_equal(lower<=.002,c['compatible']);assert np.allclose(lower,c['paired_difference_lower'],atol=1e-12)
        profiles_checked+=1
        ev=read(p)
        for row in ev.get('rows',[]):
            d=row['decoder'];pred=read(p.parent/f'prediction-{d}.json');base=pred['incumbent'];paths=sorted(p.parent.glob(f'catalog-{d}-*.json'));cats={}
            negative=[] if ev['source']=='tissue' else (list(range(998,1024)) if ev['source']=='culture' else None)
            if negative is None:
                meta=read(ROOT/'outputs/barseq3-native-repair-oct05/codebook-classes.json');names=list(dict.fromkeys(meta['names']));negative=sorted(set(names.index(n) for n,blank in zip(meta['names'],meta['blank']) if blank))
            for path in paths:
                record=read(path);cat=record['catalog'];assert len(cat['xyz'])==len(cat['words']);cat['negative_ids']=negative;cats[record['candidate']]=cat;catalogs_checked+=1
            assert len(cats)==len(prof['maps_zyx']);risk=pred['risk']
            if not risk['abstained']:
                indices=np.flatnonzero(c['compatible']);values=[comparison(cats[base],cats[int(i)]) for i in indices];identity=max(x for x,y in values);tv=max(y for x,y in values);assert abs(identity-risk['identity_risk'])<1e-12 and abs(tv-risk['composition_risk'])<1e-12;assert risk['flagged']==(identity>=.05 or tv>.05)
            assert pred['freeze_sha256']==sha(HERE/'freeze.json');prediction_checks+=1
        verified.append(dict(file=str(p.relative_to(ROOT)),sha256=sha(p)))
    out=dict(passed=True,frozen_files_verified=len(freeze['files']),case_profiles_verified=profiles_checked,prediction_records_verified=prediction_checks,catalog_records_verified=catalogs_checked,common_support_profile_bootstrap_reproduced=True,global_maximum_cardinality_matching_verified=True,composition_excludes_absent_codes=True,independent_RNA_truth=False,independent_preparation_validation=False,records=verified);(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
if __name__=='__main__':main()
