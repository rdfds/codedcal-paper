"""Retained fresh-study interpretation; no outcome-driven parameter changes."""
import collections,hashlib,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent;DEST=HERE/'tissue';ROOT=HERE.parent.parent
KINDS=['zero','common','color','outclass'];DECODERS=['source_spot','correlation','regression']
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1024**2):h.update(b)
 return h.hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def main():
 freeze=read(DEST/'freeze.json')
 for rel,digest in freeze['files'].items():assert sha(ROOT/rel)==digest,rel
 native=read(DEST/'native-state.json');assert native['completed']==15
 report=dict(native_eligibility_stage=native['stage'],original_whole_interface_gates_not_passed=True,negative_proxy_unavailable=True,native_RNA_truth_unavailable=True,full_comparative_workflow_frozen_before_pixels=True,source='distinct published mouse-brain tissue experiment; one field, five windows',freeze_sha256=sha(DEST/'freeze.json'))
 if native['stage']!='eligible':
  report.update(stage='native_eligibility_failure',cohort_comparison_supported=False,groups=native['groups']);write(DEST/'interpretation.json',report);print(json.dumps(report,indent=2));return
 state=read(DEST/'comparison-state.json');assert state['stage']=='complete' and state['completed']+state['failed']==96
 rows=[];verified=[];cases=[]
 for kind in KINDS:
  geom=read(DEST/f'{kind}-geometry.json');decisions=[]
  for dec in DECODERS:
   reference=read(DEST/'native'/f'4__{dec}.json')['catalog'];ca=collections.Counter(reference['words']);arms=[]
   for p in sorted((DEST/'outcomes').glob(f'{kind}__*__{dec}.json')):
    r=read(p);rows.append(r);verified.append(dict(file=str(p.relative_to(ROOT)),sha256=sha(p),status=r['status']))
    if r['status']!='completed':continue
    pred=r['prediction'];assert pred['negative'] is None and pred['coding']==len(pred['words'])==len(pred['xyz']);assert pred['coding']+pred['abstained']==pred['candidate_objects']
    cb=collections.Counter(pred['words']);tv=sum(abs(ca[w]/sum(ca.values())-cb[w]/sum(cb.values())) for w in ca.keys()|cb.keys())/2 if cb else None
    for m in r['metrics'].values():
     assert m['reference_calls']==len(reference['words']) and m['predicted_calls']==len(pred['words'])
     if tv is not None:assert abs(m['composition_tv']-tv)<1e-12
     assert abs(m['identity_f1']-2*m.get('correct_labels',0)/(len(reference['words'])+len(pred['words'])))<1e-12
    if r['policy']=='oracle':assert pred['xyz']==reference['xyz'] and pred['words']==reference['words']
    arms.append(r)
   conventional=[a for a in arms if a['policy'] not in ['none','oracle','global_group']]
   best=max(conventional,key=lambda a:a['metrics']['3']['identity_f1']);selected=next((a for a in arms if a['policy']==geom['selected']),None);global_arm=next((a for a in arms if a['policy']=='global_group'),None)
   sf=selected['metrics']['3']['identity_f1'] if selected else None;gf=global_arm['metrics']['3']['identity_f1'] if global_arm else None;bf=best['metrics']['3']['identity_f1']
   decisions.append(dict(decoder=dec,selected=geom['selected'],selected_f1=sf,best_conventional_policy_by_F1_evaluation_only=best['policy'],best_conventional_f1=bf,global_f1=gf,global_gain_over_best_conventional=gf-bf if gf is not None else None,selection_regret=bf-sf if sf is not None else None,rows=[dict(policy=a['policy'],f1=a['metrics']['3']['identity_f1'],composition_tv=a['metrics']['3']['composition_tv'],coding=a['prediction']['coding'],abstained=a['prediction']['abstained'],unsupported=a['prediction']['unsupported_candidates'],geometry_s=a['policy_geometry_s'],decode_s=a['decode_s']) for a in arms]))
  cases.append(dict(kind=kind,decision_cases=decisions,coding_losses=geom['coding_losses']))
 primary=[d for c in cases if c['kind']!='outclass' for d in c['decision_cases']];outclass=next(c for c in cases if c['kind']=='outclass')['decision_cases']
 full=state['failed']==0
 avoids_harm=all(d['selection_regret'] is not None and d['selection_regret']<=.02+1e-12 for d in primary)
 benefit_exists=any(d['global_gain_over_best_conventional'] is not None and d['global_gain_over_best_conventional']>=.05 for d in primary)
 catches_benefit=all(d['selected']=='global_group' for d in primary if d['global_gain_over_best_conventional'] is not None and d['global_gain_over_best_conventional']>=.05)
 diagnostic_harm=any(d['selection_regret'] is None or d['selection_regret']>.02+1e-12 for d in outclass)
 report.update(stage='interpreted',planned=96,completed=state['completed'],failed=state['failed'],cohort_comparison_supported=full,primary_rule_avoids_harm=avoids_harm,primary_distinct_benefit_exists=benefit_exists,primary_rule_catches_benefit=catches_benefit,outclass_selection_harm=diagnostic_harm,conditional_decision_confirmation=full and avoids_harm and catches_benefit,distinct_benefit_confirmation=full and avoids_harm and catches_benefit and benefit_exists,general_recommendation_established=False,cases=cases,retained_records=verified)
 write(DEST/'interpretation.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ['cases','retained_records']},indent=2))
if __name__=='__main__':main()
