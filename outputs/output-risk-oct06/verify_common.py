"""Independent common-domain bootstrap, matching and evaluation recount."""
import hashlib,json
from pathlib import Path
import numpy as np
from verify import comparison
ROOT=Path(__file__).resolve().parent.parent.parent;HERE=ROOT/'outputs/output-risk-oct06/common-domain'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    freeze=read(HERE/'freeze.json')
    for rel,h in freeze['files'].items():assert sha(ROOT/rel)==h,rel
    records=[];nc=0;npred=0;neval=0
    for path in sorted((HERE/'cases').glob('*/profile.json')):
        p=read(path);ev=read(path.parent/'evaluation.json');r=np.array(p['residual']);e=np.array(p['energy']);n=np.array(p['block_counts']);maps=np.array(p['maps_zyx']);assert np.all(maps[:,0]==0)
        use=(n>0)&(e.max(axis=0)>0);r=r[:,use];e=e[:,use];size=int(use.sum());score=r.sum(axis=1)/np.maximum(e.sum(axis=1),1e-20)
        if p['compatibility']['abstained']:assert size<12 or p['common_coverage']<.8;records.append(dict(case=path.parent.name,abstained=True));continue
        rng=np.random.default_rng(61021);indices=rng.integers(0,size,(500,size));boot=np.stack([r[:,i].sum(axis=1)/np.maximum(e[:,i].sum(axis=1),1e-20) for i in indices],axis=1);best=int(score.argmin());lower=np.quantile(boot-boot[best],.025,axis=1);compatible=lower<=.002;ties=boot<=boot.min(axis=0)+1e-6;weights=(ties/ties.sum(axis=0)).mean(axis=1)
        assert np.allclose(score,p['compatibility']['losses'],atol=1e-12);assert np.array_equal(compatible,p['compatibility']['compatible']);assert np.allclose(weights,p['selection']['weights'],atol=1e-12)
        for row in ev['rows']:
            dec=row['decoder'];pred=read(path.parent/f'prediction-{dec}.json');inc=pred['incumbent'];cats=[read(path.parent/f'catalog-{dec}-{i}.json')['catalog'] for i in range(len(maps))];nc+=len(cats)
            values=np.array([comparison(cats[inc],c) for c in cats]);maximum=values[compatible].max(axis=0);weighted=weights@values
            for key,v in [('risk',maximum),('weighted_risk',weighted)]:
                q=pred[key]
                if not q['abstained']:assert np.allclose(v,[q['identity_risk'],q['composition_risk']],atol=1e-12);assert q['flagged']==(v[0]>=.05 or v[1]>.05)
                else:assert len([w for w in cats[inc]['words'] if w not in cats[inc]['negative_ids']])<10
                npred+=1
            ref=read(path.parent/f'reference-{dec}.json')['catalog'];actual=comparison(ref,cats[inc]);a=row['actual_consistency_error'];assert np.allclose(actual,[a['identity_loss'],a['composition_tv']],atol=1e-12);assert row['failed_consistency']==(actual[0]>.05 or actual[1]>.05);neval+=1
        records.append(dict(case=path.parent.name,sha256=sha(path.parent/'evaluation.json'),common_coverage=p['common_coverage']))
    assert len(records)==13,len(records)
    out=dict(passed=True,frozen_files=len(freeze['files']),profiles=len(records),candidate_catalogs=nc,forecast_statistic_records=npred,evaluation_records=neval,bootstrap_and_weighted_ties_reproduced=True,canonical_decoder_canvas_verified=True,maximum_cardinality_matching_verified=True,native_reference_recounted=True,independent_RNA_truth=False,independent_preparation_validation=False,records=records);(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
if __name__=='__main__':main()
