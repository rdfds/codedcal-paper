"""Independent follow-up recount and global tissue assignment sensitivity."""
import collections,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    frozen=[]
    for spec in ['freeze.json','local-trace/freeze.json','tissue/freeze.json','fullfield-tissue/comparison-freeze.json']:
        for rel,digest in read(HERE/spec)['files'].items():assert sha(ROOT/rel)==digest,rel;frozen.append(rel)
    core=ROOT/'outputs/barseq3-native-repair-oct05/global_repair.py';assert sha(core)=='82407386fccc55da8b745fd5a31c83ed27d41e4c4fcb2eb1b3267d5c14df4fdc'
    records=[];sens=[]
    paths=list((HERE/'fullfield-tissue/outcomes').glob('*.json'))
    assert len(paths)==96
    if (HERE/'volumetric-model/state.json').exists():paths+=list((HERE/'volumetric-model').glob('color__*.json'))
    for p in paths:
        r=read(p);assert r['status']=='completed';pred=r['prediction'];ref=read(HERE/'fullfield-tissue/native'/f"4__{r['decoder']}.json")['catalog'];a=np.array(ref['xyz']).reshape(-1,3);b=np.array(pred['xyz']).reshape(-1,3);aw=np.array(ref['words']);bw=np.array(pred['words']);assert len(b)==len(bw)==pred['coding'] and pred['negative'] is None;assert pred['coding']+pred['abstained']==pred['candidate_objects'];ca=collections.Counter(aw);cb=collections.Counter(bw);tv=sum(abs(ca[w]/len(a)-cb[w]/len(b)) for w in ca.keys()|cb.keys())/2 if len(b) else None
        dist=cdist(a,b)
        for radius in [1.5,3,5]:
            primary=r['metrics'][str(radius)];assert primary['reference_calls']==len(a) and primary['predicted_calls']==len(b);assert abs(primary['identity_f1']-2*primary['correct_labels']/(len(a)+len(b)))<1e-12
            if tv is not None:assert abs(primary['composition_tv']-tv)<1e-12
            valid=dist<=radius;maximum=int((maximum_bipartite_matching(csr_matrix(valid.astype(int)),perm_type='column')>=0).sum());reward=(min(len(a),len(b))+1)*radius;cost=np.zeros((len(a)+len(b),len(a)+len(b)));cost[:len(a),:len(b)]=np.where(valid,dist-reward,reward);i,j=linear_sum_assignment(cost);ok=(i<len(a))&(j<len(b));i,j=i[ok],j[ok];ok=valid[i,j];i,j=i[ok],j[ok];assert len(i)==maximum;f1=2*int((aw[i]==bw[j]).sum())/(len(a)+len(b));sens.append(dict(file=str(p.relative_to(ROOT)),kind=r['kind'],policy=r['policy'],decoder=r['decoder'],radius=radius,primary_f1=primary['identity_f1'],global_f1=f1,difference=f1-primary['identity_f1']))
        if r['metrics']['3']['identity_f1']==1:assert pred['xyz']==ref['xyz'] and pred['words']==ref['words']
        records.append(dict(file=str(p.relative_to(ROOT)),sha256=sha(p)))
    for folder in ['object-score','source-anchor-score']:
        assert read(HERE/folder/'state.json')['stage']=='complete'
        for kind in ['zero','common','color','outclass']:
            p=HERE/folder/f'{kind}.json';r=read(p);records.append(dict(file=str(p.relative_to(ROOT)),sha256=sha(p)))
    report=dict(passed=True,tissue_comparison_records=len(paths),decision_cases=8,records=records,frozen_files_verified=len(set(frozen)),unchanged_core_sha256=sha(core),identity_counts_composition_verified=True,matching_radius_comparisons=len(sens),matching_maximum_observed_difference=max(abs(s['difference']) for s in sens),matching_different_cases=sum(abs(s['difference'])>1e-12 for s in sens),matching_uncertainty_not_exhaustive=True,independent_RNA_truth=False,independent_complete_workflow_confirmation=False,matching_sensitivity=sens)
    (HERE/'followup-verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['records','matching_sensitivity']},indent=2))
if __name__=='__main__':main()
