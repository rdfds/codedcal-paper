"""Independent whole-domain matching check for both external cohorts."""
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from condition import HERE,read,write,sha
from verify_results import spatial_match

def main():
    classes=read(HERE/'codebook-classes.json');names=list(dict.fromkeys(classes['names']))
    label=np.array([names.index(n) for n in classes['names']]);rows=[]
    for cohort in ['recovery','fit-support']:
        summary=read(HERE/cohort/'summary.json')
        assert summary['counts'].get('missing',0)==0 and summary['counts'].get('not_run',0)==0
        for index in [0,1]:
            for p in sorted((HERE/cohort/f'field-{index}').glob('*__*.json')):
                arm=read(p)
                if arm['status']!='completed':continue
                ref=read(HERE/'recovery'/f'field-{index}'/f"reference-{arm['decoder']}.json")
                a=np.array(ref['xy'],float).reshape(-1,2);b=np.array(arm['xy'],float).reshape(-1,2)
                for radius in [1.5,3,5]:
                    dist=cdist(a,b);edges=dist<=radius;n,m=dist.shape
                    cardinality=int((maximum_bipartite_matching(csr_matrix(edges),perm_type='column')>=0).sum()) if n and m else 0
                    pairs,diag=spatial_match(a,b,radius)
                    assert diag['skipped_components']==0 and len(pairs)==cardinality
                    # Whole-domain rectangular assignment with one unmatched dummy per reference.
                    reward=(min(n,m)+1)*radius
                    cost=np.zeros((n,m+n));cost[:,:m]=np.where(edges,dist-reward,reward)
                    rr,cc=linear_sum_assignment(cost)
                    valid=(cc<m)
                    rr,cc=rr[valid],cc[valid]
                    valid=edges[rr,cc];rr,cc=rr[valid],cc[valid]
                    assert len(rr)==cardinality
                    local_distance=float(np.linalg.norm(a[pairs[:,0]]-b[pairs[:,1]],axis=1).sum())
                    global_distance=float(dist[rr,cc].sum())
                    assert abs(local_distance-global_distance)<1e-8
                    correct=int((label[np.array(ref['words'],int)[rr]]==label[np.array(arm['words'],int)[cc]]).sum())
                    primary=arm['metrics'][str(radius)]['correct_gene_classes']
                    rows.append(dict(cohort=cohort,index=index,kind=arm['kind'],policy=arm['policy'],decoder=arm['decoder'],radius=radius,
                        spatial_matches=cardinality,global_identity_correct=correct,primary_identity_correct=primary,
                        tied_matching_identity_difference=correct-primary,total_distance=global_distance,
                        global_identity_f1=2*correct/(n+m) if n+m else 0.,primary_identity_f1=arm['metrics'][str(radius)]['gene_class_identity_f1']))
    write(HERE/'assignment-verification.json',dict(rows=rows,completed_arms=len(rows)//3,
        global_cardinality_and_distance_verified=True,
        equal_cost_identity_differences=sum(r['tied_matching_identity_difference']!=0 for r in rows),
        largest_observed_equal_cost_f1_change=max(abs(r['global_identity_f1']-r['primary_identity_f1']) for r in rows),
        script_sha256=sha(__file__),scope='Independent sparse cardinality and whole-domain distance assignment; tied identity differences retained, not RNA truth'))
    print('Verified',len(rows)//3,'arms; equal-cost identity differences:',sum(r['tied_matching_identity_difference']!=0 for r in rows),flush=True)

if __name__=='__main__':main()
