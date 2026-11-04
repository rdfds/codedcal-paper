"""Evaluate model adequacy after saved forecasts; never repair by truth selection."""
import json
from pathlib import Path
import numpy as np
from verify import comparison
ROOT=Path(__file__).resolve().parent.parent.parent;HERE=ROOT/'outputs/output-risk-oct06/common-domain'
def read(p):return json.loads(p.read_text())
def main():
 records=[]
 for path in sorted((HERE/'cases').glob('*/evaluation.json')):
  ev=read(path)
  if ev['stage']!='complete':continue
  prof=read(path.parent/'profile.json');loss=np.array(prof['compatibility']['losses']);weights=np.array(prof['selection']['weights']);labels=prof['labels']
  for row in ev['rows']:
   dec=row['decoder'];ref=read(path.parent/f'reference-{dec}.json')['catalog'];cats=[read(path.parent/f'catalog-{dec}-{k}.json')['catalog'] for k in range(len(loss))];error=np.array([comparison(ref,c) for c in cats]);bad=(error[:,0]>.05)|(error[:,1]>.05);min_error=np.min(np.maximum(error[:,0],error[:,1]));acceptable=np.flatnonzero(~bad);best=int(loss.argmin());inc=row['incumbent'];order=np.argsort(loss)
   records.append(dict(case=ev['case'],decoder=dec,incumbent_labels=labels[inc],incumbent_loss=float(loss[inc]),minimum_catalog_error_in_bank=float(min_error),best_fit_labels=labels[best],best_fit_native_error=error[best].tolist(),best_native_consistent_loss=None if not len(acceptable) else float(loss[acceptable].min()),strict_fit_preference_for_native_inconsistent=bool(len(acceptable) and bad[best] and loss[best]<loss[acceptable].min()-1e-6),bootstrap_weight_on_native_inconsistent=float(weights[bad].sum()),incumbent_native_error=error[inc].tolist(),highest_weight_labels=labels[int(weights.argmax())],top_five_fit_candidates=[dict(labels=labels[k],loss=float(loss[k]),native_error=error[k].tolist(),weight=float(weights[k])) for k in order[:5]],evaluation_only=True))
 out=dict(stage='model_adequacy_diagnosis',cases=len(set(r['case'] for r in records)),records=records,truth_used_for_diagnosis_only=True,no_outcome_selected_forecast_repair=True,independent_RNA_truth=False);(HERE/'diagnosis.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps([ {k:v for k,v in r.items() if k!='top_five_fit_candidates'} for r in records],indent=2))
if __name__=='__main__':main()
