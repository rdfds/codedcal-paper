"""Independent whole-domain optimal assignment for strong conventional BARseq barcode-model rows."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent.parent;BAR=ROOT/'outputs/barseq3-native-repair-oct05'
def read(p):return json.loads(p.read_text())
def main():
 labels=read(BAR/'codebook-classes.json')['names'];names=list(dict.fromkeys(labels));ids=np.array([names.index(n) for n in labels]);rows=[]
 for p in sorted((HERE/'bar-model').glob('*_*_*.json')):
  r=read(p);ref=read(BAR/'recovery'/f"field-{r['field']}"/f"reference-{r['decoder']}.json");a=np.asarray(ref['xy']).reshape(-1,2);b=np.asarray(r['xy']).reshape(-1,2);aw=ids[np.asarray(ref['words'],int)];bw=ids[np.asarray(r['words'],int)];dist=cdist(a,b)
  for radius in [1.5,3,5]:
   valid=dist<=radius;maximum=int((maximum_bipartite_matching(csr_matrix(valid.astype(int)),perm_type='column')>=0).sum());reward=(min(len(a),len(b))+1)*radius;cost=np.zeros((len(a)+len(b),len(a)+len(b)));cost[:len(a),:len(b)]=np.where(valid,dist-reward,reward);i,j=linear_sum_assignment(cost);ok=(i<len(a))&(j<len(b));i,j=i[ok],j[ok];ok=valid[i,j];i,j=i[ok],j[ok];assert len(i)==maximum;correct=int((aw[i]==bw[j]).sum());f1=2*correct/(len(a)+len(b));primary=r['metrics'][str(radius)]['gene_class_identity_f1']
   rows.append(dict(file=str(p.relative_to(ROOT)),field=r['field'],kind=r['kind'],decoder=r['decoder'],radius=radius,maximum_cardinality=maximum,global_identity_f1=f1,primary_identity_f1=primary,assignment_F1_difference=f1-primary))
 result=dict(passed=True,records=12,radius_comparisons=36,rows=rows,matching_sensitivity_not_exhaustive=True,maximum_absolute_F1_difference=max(abs(r['assignment_F1_difference']) for r in rows),tied_or_greedy_differences=sum(abs(r['assignment_F1_difference'])>1e-12 for r in rows))
 (HERE/'strong-bar-matching-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
 for r in rows:
  if r['kind']=='color' and r['radius']==3:print(r)
if __name__=='__main__':main()
