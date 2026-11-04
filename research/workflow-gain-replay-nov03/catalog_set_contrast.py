"""Post-hoc set comparison; raw sorted row order is not spot identity.

Preserve the original ordered-row diagnostic. Tiny centroid rounding can reorder
spots when a lexicographic sort changes a tied ordinate. Use the existing exact
cardinality matcher at the frozen engineering tolerance, separately by class.
No image read, decoder prediction, fitting, threshold change or rescue re-score.
"""
import fcntl
from replay_common import *
from codedcal.workflow_validation import spatial_match


def main():
    f=load(H/'freeze.json');out=H/'run'
    lock=(out/'collector.lock').open('r');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert load(out/'state.json')['stage']=='completed'
    original=load(H/'posthoc-contrast.json');rows=[]
    for e in f['cases']:
        a=load(out/'historical'/(e['id']+'.json'));b=load(out/'anchored'/(e['id']+'.json'))
        ac,bc=Catalog.from_dict(a['catalog']),Catalog.from_dict(b['catalog'])
        matched=0;unresolved=0;maximum=0.
        for word in np.union1d(ac.classes,bc.classes):
            aa=ac.points[ac.classes==word];bb=bc.points[bc.classes==word]
            pairs,info=spatial_match(aa,bb,f['coordinate_engineering_tolerance'])
            matched+=len(pairs);unresolved+=info['unresolved_components']
            if len(pairs):maximum=max(maximum,float(np.max(np.linalg.norm(aa[pairs[:,0]]-bb[pairs[:,1]],axis=1))))
        equivalent=bool(unresolved==0 and matched==len(ac.points)==len(bc.points))
        old=next(r for r in original['rows'] if r['case_id']==e['id'])
        rows.append(dict(case_id=e['id'],native=e['challenge'] is None,matched_spots=matched,
            historical_spots=len(ac.points),anchored_spots=len(bc.points),unresolved_components=unresolved,
            same_class_catalog_sets_equivalent=equivalent,maximum_matched_distance_pixels=maximum,
            raw_rowwise_equivalent=bool(old['catalog_counts_and_ordered_classes_identical'] and old['maximum_coordinate_difference_pixels']<=1e-12)))
    result=dict(post_hoc_diagnostic=True,rows=rows,cases=26,integer_cases=24,
        catalog_sets_equivalent=sum(r['same_class_catalog_sets_equivalent'] for r in rows),
        equivalent_sets_with_reordered_rows=sum(r['same_class_catalog_sets_equivalent'] and not r['raw_rowwise_equivalent'] for r in rows),
        maximum_matched_distance_pixels=max(r['maximum_matched_distance_pixels'] for r in rows),
        engineering_tolerance=f['coordinate_engineering_tolerance'],
        molecular_scoring_unchanged=True,preserved_original_ordered_diagnostic='posthoc-contrast.json',
        new_decoder_executions=0,new_optimizations=0)
    target=H/'catalog-set-contrast.json'
    if target.exists():assert load(target)==result
    else:write_json(target,result)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))


if __name__=='__main__':main()
