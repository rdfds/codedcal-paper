"""Post-hoc diagnostic of completed replay catalogs; no images or predictions.

The rescue endpoint is already frozen/adjudicated. This separate diagnostic
checks whether gains actually changed while catalogs and scores stayed stable.
"""
import fcntl
from replay_common import *


def main():
    f=load(H/'freeze.json');out=H/'run'
    lock=(out/'collector.lock').open('r');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert load(out/'state.json')['stage']=='completed'
    rows=[]
    for e in f['cases']:
        a=load(out/'historical'/(e['id']+'.json'));b=load(out/'anchored'/(e['id']+'.json'))
        assert a['stage']==b['stage']=='completed'
        ac,bc=Catalog.from_dict(a['catalog']),Catalog.from_dict(b['catalog'])
        labels_equal=np.array_equal(ac.classes,bc.classes) and ac.points.shape==bc.points.shape
        error=float(np.max(np.abs(ac.points-bc.points))) if labels_equal and len(ac.points) else (0. if labels_equal else None)
        aa,bb=np.asarray(a['gains_used']),np.asarray(b['gains_used'])
        offsets=np.log(aa)-np.log(bb)
        rows.append(dict(case_id=e['id'],policy=e['policy'],native=e['challenge'] is None,
            catalog_counts_and_ordered_classes_identical=labels_equal,maximum_coordinate_difference_pixels=error,
            all_assessment_fields_identical=a['assessment']==b['assessment'],quality_identical=a['quality']==b['quality'],
            gains_changed=not np.array_equal(aa,bb),common_gain_log_multiplier=b['gain_gauge']['log_common_gain_multiplier'],
            centered_log_gain_difference_range=float(np.ptp(offsets))))
    result=dict(post_hoc_diagnostic=True,rows=rows,cases=26,integer_cases=24,
        catalogs_equivalent_at_frozen_coordinate_tolerance=sum(r['catalog_counts_and_ordered_classes_identical'] and r['maximum_coordinate_difference_pixels']<=1e-12 for r in rows),
        maximum_coordinate_difference_pixels=max(r['maximum_coordinate_difference_pixels'] for r in rows if r['maximum_coordinate_difference_pixels'] is not None),
        identical_assessments=sum(r['all_assessment_fields_identical'] for r in rows),
        changed_gain_vectors=sum(r['gains_changed'] for r in rows),
        new_decoder_executions=0,new_optimizations=0)
    target=H/'posthoc-contrast.json'
    if target.exists():assert load(target)==result
    else:write_json(target,result)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))


if __name__=='__main__':main()
