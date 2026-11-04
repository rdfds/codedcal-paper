"""Separate weighted forecasts; neither truth nor known offset enters prediction."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent.parent.parent;HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal import output_risk_weighted as W

def read(p):return json.loads(p.read_text())
def write(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def main():
    f=read(HERE/'weighted-freeze.json')
    for rel,digest in f['files'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
    for path in sorted((HERE/'cases').glob('*/profile.json')):
        prof=read(path);profile=dict(residual=np.asarray(prof['residual']),energy=np.asarray(prof['energy']),block_counts=np.asarray(prof['block_counts']),coverage=prof['common_coverage']);selection=W.selection_weights(profile);dest=HERE/'weighted/cases'/path.parent.name;updates={}
        for pp in sorted(path.parent.glob('prediction-*.json')):
            pred=read(pp);dec=pred['decoder'];cats=[]
            for k in range(len(prof['maps_zyx'])):
                cat=read(path.parent/f'catalog-{dec}-{k}.json')['catalog']
                if path.parent.name.startswith('tissue'):negative=[]
                elif path.parent.name.startswith('culture'):negative=list(range(998,1024))
                else:
                    meta=read(ROOT/'outputs/barseq3-native-repair-oct05/codebook-classes.json');names=list(dict.fromkeys(meta['names']));negative=sorted(set(names.index(n) for n,b in zip(meta['names'],meta['blank']) if b))
                cat['negative_ids']=negative;cats.append(cat)
            risk=W.weighted_output_risk(cats,selection,pred['incumbent'],np.asarray(prof['maps_zyx']));old=pred['risk']
            if not risk['abstained']:
                for key in ['ordinary_multistart_identity_disagreement','ordinary_multistart_composition_disagreement','unfiltered_bank_identity_disagreement','unfiltered_bank_composition_disagreement']:risk[key]=old[key]
            write(dest/f'prediction-{dec}.json',dict(decoder=dec,risk=risk,selection=selection,development=True,posthoc_statistic_repair=True,posterior_RNA_error_probability=False));updates[dec]=risk
        # Truth is consulted only after separate forecasts have been saved.
        original=read(path.parent/'evaluation.json');rows=[r|dict(risk=updates[r['decoder']]) for r in original['rows']];write(dest/'evaluation.json',original|dict(rows=rows,posthoc_statistic_repair=True));print('WEIGHTED',path.parent.name,{d:r.get('identity_risk') for d,r in updates.items()},flush=True)
if __name__=='__main__':main()
