"""Recount original-domain weighted forecasts without calling the risk library."""
import hashlib,json
from pathlib import Path
import numpy as np
from verify import comparison
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent

def read(p):return json.loads(p.read_text())
def main():
 f=read(HERE/'weighted-freeze.json')
 for rel,h in f['files'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==h
 n=0
 for p in sorted((HERE/'cases').glob('*/profile.json')):
  prof=read(p);r=np.asarray(prof['residual']);e=np.asarray(prof['energy']);use=(np.array(prof['block_counts'])>0)&(e.max(axis=0)>0);r=r[:,use];e=e[:,use];size=int(use.sum());rng=np.random.default_rng(61021);indices=rng.integers(0,size,(500,size));loss=np.stack([r[:,i].sum(axis=1)/np.maximum(e[:,i].sum(axis=1),1e-20) for i in indices],axis=1);ties=loss<=loss.min(axis=0)+1e-6;weights=(ties/ties.sum(axis=0)).mean(axis=1)
  for pp in sorted(p.parent.glob('prediction-*.json')):
   original=read(pp);dec=original['decoder'];new=read(HERE/'weighted/cases'/p.parent.name/f'prediction-{dec}.json');assert np.allclose(weights,new['selection']['weights'],atol=1e-12);cats=[read(p.parent/f'catalog-{dec}-{k}.json')['catalog'] for k in range(len(weights))]
   if p.parent.name.startswith('culture'):negative=list(range(998,1024))
   elif p.parent.name.startswith('tissue'):negative=[]
   else:
    meta=read(ROOT/'outputs/barseq3-native-repair-oct05/codebook-classes.json');names=list(dict.fromkeys(meta['names']));negative=sorted(set(names.index(n) for n,b in zip(meta['names'],meta['blank']) if b))
   for cat in cats:cat['negative_ids']=negative
   inc=original['incumbent'];coding=sum(w not in negative for w in cats[inc]['words']);q=new['risk']
   if coding<10:assert q['abstained']
   else:
    vals=np.array([comparison(cats[inc],c) for c in cats]);v=weights@vals;assert np.allclose(v,[q['identity_risk'],q['composition_risk']],atol=1e-12);assert q['flagged']==(v[0]>=.05 or v[1]>.05)
   n+=1
 assert n==20;out=dict(passed=True,weighted_forecasts_recounted=n,independent_tie_weight_and_global_assignment_recount=True,independent_RNA_truth=False);(HERE/'weighted/verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
