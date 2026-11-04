"""Artificial/pre-outcome controls; no geometric preparation outcomes inspected."""
import json
import numpy as np
from methods import estimate, interval, positive_parts, product_root, contains, geometric_target, arithmetic_target
from run import HERE, write


def main():
    checks={}
    a=.9; f=.7
    baseline=np.array([[a,1.],[a,1.]])
    changed=np.array([[f*a,1.],[f*a,f*f]])
    np.testing.assert_allclose(geometric_target(changed),a,atol=1e-15)
    assert abs(arithmetic_target(changed)-a)>.01
    groups=np.broadcast_to(changed,(32,2,2)).copy()
    result=estimate(groups,'geometric')
    assert result['usable'] and contains(result['intervals'],a)
    checks['uniform_kernel_factor_cancels']=True
    wA,wB=np.array([3.,8.]),np.array([4.,4.]);factor=np.array([.5,2.])
    variable=np.array([[np.sum(wA*factor*a),wA.sum()],[np.sum(wB*factor*a),np.sum(wB*factor**2)]])
    assert abs(geometric_target(variable)-a)>.01
    checks['source_dependent_kernel_counterexample']=True
    weights=np.array([[1.,3.],[3.,1.]])
    source=np.array([.7,.95]);mixed=np.column_stack([weights@source,weights.sum(1)])
    assert abs(geometric_target(mixed)-arithmetic_target(mixed))>1e-4
    checks['heterogeneous_window_weights_change_target']=True
    negatives=np.broadcast_to(np.array([[-1.,1.],[-1.,1.]]),(16,2,2)).copy()
    rejected=estimate(negatives,'geometric')
    assert not rejected['usable'] and rejected['ratio'] is None and rejected['intervals']==[]
    checks['two_negative_ratios_cannot_earn_call']=True
    badref=np.broadcast_to(np.array([[-1.,-1.],[-1.,-1.]]),(16,2,2)).copy()
    rejected=estimate(badref,'geometric')
    assert rejected['ratio'] is None and not rejected['usable']
    checks['negative_references_rejected']=True
    assert product_root([interval(1,4,True,True)],[interval(4,9,True,True)])==[interval(2.,6.,True,True)]
    assert product_root([interval(0,None)],[interval(1,2,True,True)])==[interval(0.,None)]
    assert product_root([], [interval(1,2)])==[]
    assert product_root(None,[]) is None
    checks['product_set_endpoints']=True
    disjoint=[interval(None,1.,False,True),interval(2.,None,True,False)]
    positive=positive_parts(disjoint)
    assert product_root(positive,positive)==[interval(0.,None)]
    assert not contains(product_root(positive,positive),0)
    checks['disjoint_sets_union_without_finite_rescue']=True
    for groups in [np.zeros((8,2,2)),negatives,badref,np.broadcast_to(baseline,(16,2,2))]:
        for variant in ('arithmetic','geometric'):
            json.dumps(estimate(groups,variant),allow_nan=False)
    checks['strict_json_degenerate_cases']=True
    write(HERE/'controls.json',dict(all_pass=True,checks=checks,pre_outcome=True))
    print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
