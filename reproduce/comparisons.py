"""Recount complete cohorts and negative comparisons without fitting or decoding."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
def read(name):return json.loads((ROOT/name).read_text())

def main():
 d=read('research/workflow-qualification-nov01/adjudicated-summary.json')
 assert d['complete'] and len(d['rows'])==d['original_planned_workflows']==32
 truth=[r for r in d['rows'] if r['truth']['decision']!='unsupported']
 failed=[r for r in truth if r['truth']['decision']=='reject']
 adequate=[r for r in truth if r['truth']['decision']=='accept']
 assert (len(truth),len(failed),len(adequate))==(10,9,1)
 for name in ('mixed_conventional','readout_aware'):
  c=d['confusion'][name]
  for field,cohort,decision in [('failed_rejected',failed,'reject'),('failed_accepted',failed,'accept'),('adequate_rejected',adequate,'reject'),('adequate_accepted',adequate,'accept')]:
   assert c[field]==sum(r[name]['decision']==decision for r in cohort)
 assert not d['empirical_utility_comparison_supported']
 r=read('research/workflow-validation-completion-nov03/summary.json')
 assert r['complete'] and r['complete_integer_matrix'] and len(r['rows'])==4
 eligible=[x for x in r['rows'] if x['native_eligible']]
 assert len(eligible)==r['supported_heldout_cases']==2
 assert sum(x['truth']['decision']=='accept' for x in eligible)==r['adequate_cases']==0
 assert all(x['truth']['decision']=='reject' and x['mixed_conventional']['decision']=='reject' and x['mechanism_weighted']['decision']=='reject' for x in eligible)
 assert not r['general_utility_supported'] and not r['comparison_rates_estimable']
 g=read('research/workflow-gain-replay-nov03/summary.json')
 c=read('research/workflow-gain-replay-nov03/catalog-set-contrast.json')
 assert g['complete'] and g['completed_predictions']==g['actual_decoder_predictions_started']==52
 assert len(g['case_changes'])==24
 assert sum(x['before_failed'] and not x['after_failed'] for x in g['case_changes'])==g['failures_to_passes']==0
 assert sum(not x['before_failed'] and x['after_failed'] for x in g['case_changes'])==g['passes_to_failures']==0
 assert g['development_workflows_rescued']==0
 assert len(c['rows'])==c['catalog_sets_equivalent']==26
 assert all(x['same_class_catalog_sets_equivalent'] for x in c['rows'])
 assert sum(not x['raw_rowwise_equivalent'] for x in c['rows'])==c['equivalent_sets_with_reordered_rows']==25
 assert max(x['maximum_matched_distance_pixels'] for x in c['rows'])==c['maximum_matched_distance_pixels']<=c['engineering_tolerance']
 out=dict(passed=True,development_workflows=32,supported_development=10,adequate_development=1,
          complete_integer_pairs=24,assigned_policies=4,eligible_policies=2,adequate_confirmation=0,
          decoder_predictions=52,gain_case_changes=24,catalog_pairs=26,rescues=0,
          strict_superiority_supported=False,new_outcome_study=False)
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
