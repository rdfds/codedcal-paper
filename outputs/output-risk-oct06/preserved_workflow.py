"""Keep the caller's strong author workflow unchanged while diagnosing risk."""
import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent.parent.parent;HERE=ROOT/'outputs/output-risk-oct06';sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal import output_risk as O
from codedcal import output_risk_weighted as W
def read(p):return json.loads(p.read_text())
def write(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def main():
 freeze=read(HERE/'preserved-workflow-freeze.json')
 for rel,h in freeze['files'].items():assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==h,rel
 for path in sorted((HERE/'common-domain/cases').glob('*/profile.json')):
  p=read(path);source=path.parent.name.split('_')[0];kind=path.parent.name.split('_')[-1];label='author_selected' if source=='bar' else ('author_volume_selected' if source=='tissue' and kind=='color' else 'author_multistart');inc=[k for k,ls in enumerate(p['labels']) if label in ls];assert len(inc)==1,(path,label);inc=inc[0];dest=HERE/'preserved-workflow/cases'/path.parent.name;preds={};cats_by_dec={}
  for file in sorted(path.parent.glob('prediction-*.json')):
   original=read(file);dec=original['decoder'];cats=[read(path.parent/f'catalog-{dec}-{k}.json')['catalog'] for k in range(len(p['labels']))];risk=O.output_risk(cats,p['compatibility'],inc,p['conventional_indices']);weighted=W.weighted_output_risk(cats,p['selection'],inc,np.asarray(p['maps_zyx']))
   if not weighted['abstained']:
    for key in ['ordinary_multistart_identity_disagreement','ordinary_multistart_composition_disagreement','unfiltered_bank_identity_disagreement','unfiltered_bank_composition_disagreement']:weighted[key]=risk[key]
   # The saved margin describes another incumbent in some cases. Do not reuse it.
   pred=original|dict(incumbent=inc,incumbent_label=label,risk=risk,weighted_risk=weighted,baseline_fit_loss=p['compatibility']['losses'][inc],baseline_fit_gap=p['compatibility']['losses'][inc]-min(p['compatibility']['losses']),baseline_decoder_confidence=None,baseline_quality=read(path.parent/f'catalog-{dec}-{inc}.json')['quality'],confidence_comparator_not_recomputed=True)
   write(dest/f'prediction-{dec}.json',pred);preds[dec]=pred;cats_by_dec[dec]=cats
  # All decision forecasts for the case exist before evaluation-only assignments.
  ev=read(path.parent/'evaluation.json');rows=[]
  for dec,pred in preds.items():
   ref=read(path.parent/f'reference-{dec}.json')['catalog'];actual=O.catalog_distance(ref,cats_by_dec[dec][inc]);rows.append(pred|dict(reference_calls=len(ref['words']),incumbent_calls=len(cats_by_dec[dec][inc]['words']),actual_consistency_error=actual,failed_consistency=actual['identity_loss']>.05 or actual['composition_tv']>.05))
  write(dest/'evaluation.json',ev|dict(rows=rows,posthoc_preserved_workflow_sensitivity=True,decoder_confidence_comparator_not_available=True));print('PRESERVED',path.parent.name,[(r['decoder'],r['actual_consistency_error']['identity_loss']) for r in rows],flush=True)
if __name__=='__main__':main()
