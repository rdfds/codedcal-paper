"""Extended comparator assessment; does not overwrite frozen study outcomes."""
import json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent
def read(p):return json.loads(p.read_text())
def main():
    state=read(HERE/'volumetric-model/state.json');assert state['stage']=='complete' and state['completed']==3 and state['failed_starts']==0
    rows=[]
    for dec in ['source_spot','correlation','regression']:
        strong=read(HERE/'volumetric-model'/f'color__{dec}.json');group=read(HERE/'fullfield-tissue/outcomes'/f'color__global_group__{dec}.json');earlier=[read(p) for p in (HERE/'fullfield-tissue/outcomes').glob(f'color__*__{dec}.json') if 'global_group' not in p.name and 'oracle' not in p.name];best=max(earlier+[strong],key=lambda r:r['metrics']['3']['identity_f1']);f=best['metrics']['3']['identity_f1'];g=group['metrics']['3']['identity_f1'];rows.append(dict(decoder=dec,strongest_tested_policy_by_outcome_for_evaluation_only=best['policy'],strongest_tested_f1_radius3=f,group_f1_radius3=g,gain=g-f,passes_predeclared_five_point_benefit=g-f>=.05,full_volume_geometry_s=strong['policy_geometry_s'],decode_s=strong['decode_s'],composition_tv=strong['metrics']['3']['composition_tv'],coding=strong['prediction']['coding'],abstained=strong['prediction']['abstained'],unsupported=strong['prediction']['unsupported_candidates'],source_file=str((HERE/'volumetric-model'/f'color__{dec}.json').relative_to(ROOT))))
    out=dict(stage='interpreted',extended_tissue_comparator_arms=99,failed=0,development=True,original_frozen_tissue_native_eligibility_failed=True,independent_complete_workflow_confirmation=False,independent_RNA_truth=False,distinct_tissue_benefit_at_five_points=False,largest_tissue_color_gain=max(r['gain'] for r in rows),new_observable_validated=False,general_recommendation_established=False,author_and_group_have_same_physical_codebook=True,culture_strong_multistart_ties_group=True,bar_field_zero_strong_model_f1_radius3=[.9710144927536232,.9787234042553191],bar_field_one_strongest_historical_author_f1_radius3=[.029411764705882353,.008633093525179856],cases=rows)
    (HERE/'strong-comparator-summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
