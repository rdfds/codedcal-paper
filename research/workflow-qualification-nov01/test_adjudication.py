import numpy as np
from adjudicate import adjudicate_trial


def row(n,candidate_status='failed',coverage=0.):
    points=[[0,i,0] for i in range(n)]
    good=dict(status='completed',catalog=dict(points=points,classes=[0]*n),coverage=1,corrections=np.zeros((2,3)).tolist(),diagnostics={})
    a=dict(good,coverage=coverage) if candidate_status=='completed' else dict(status='failed',error='known control')
    return dict(runs=[a,good],assessment=dict(status='workflow_failed',failed=True),challenge=dict(maps=np.zeros((2,3)).tolist()))


def test_small_reference_precedes_candidate_execution_failure():
    a=adjudicate_trial(row(3),[]);assert a['assessment']['failed'] is None and a['eligibility_adjudication_changed']


def test_small_reference_precedes_output_support_failure():
    a=adjudicate_trial(row(3,'completed'),[]);assert a['assessment']['reason']=='insufficient_reference_calls'


def test_sufficient_reference_keeps_real_candidate_failure():
    a=adjudicate_trial(row(12),[]);assert a['assessment']['failed'] is True


def test_sufficient_supported_identical_catalog_is_accepted():
    a=adjudicate_trial(row(12,'completed',1),[]);assert a['assessment']['failed'] is False
