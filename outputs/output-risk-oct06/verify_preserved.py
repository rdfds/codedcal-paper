"""Recount preserved-incumbent forecasts against immutable common-domain records."""
import hashlib,json
from pathlib import Path
import numpy as np
from verify import comparison
ROOT=Path(__file__).resolve().parent.parent.parent;HERE=ROOT/'outputs/output-risk-oct06'
def read(p):return json.loads(p.read_text())
def main():
 f=read(HERE/'preserved-workflow-freeze.json')
 for rel,h in f['files'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==h
 records=[]
 for path in sorted((HERE/'preserved-workflow/cases').glob('*/evaluation.json')):
  ev=read(path);source=HERE/'common-domain/cases'/path.parent.name;p=read(source/'profile.json');weights=np.array(p['selection']['weights']);compatible=np.array(p['compatibility']['compatible'],bool)
  for r in ev['rows']:
   dec=r['decoder'];pred=read(path.parent/f'prediction-{dec}.json');inc=pred['incumbent'];assert pred['incumbent_label'] in p['labels'][inc];cats=[read(source/f'catalog-{dec}-{k}.json')['catalog'] for k in range(len(weights))];vals=np.array([comparison(cats[inc],c) for c in cats]);coding=sum(w not in cats[inc]['negative_ids'] for w in cats[inc]['words'])
   for key,v in [('risk',vals[compatible].max(axis=0)),('weighted_risk',weights@vals)]:
    risk=pred[key]
    if coding<10:assert risk['abstained']
    else:assert np.allclose([risk['identity_risk'],risk['composition_risk']],v,atol=1e-12);assert risk['flagged']==(v[0]>=.05 or v[1]>.05)
   ref=read(source/f'reference-{dec}.json')['catalog'];actual=comparison(ref,cats[inc]);assert np.allclose(actual,[r['actual_consistency_error']['identity_loss'],r['actual_consistency_error']['composition_tv']],atol=1e-12);records.append(dict(case=ev['case'],decoder=dec,incumbent=inc,coding=coding))
 assert len(records)==20
 out=dict(passed=True,decoder_records=20,forecast_statistics=40,native_assignments_recounted=True,workflow_choice_by_label_not_truth=True,confidence_comparator_not_recomputed=True,independent_RNA_truth=False,records=records);(HERE/'preserved-workflow/verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
if __name__=='__main__':main()
