import numpy as np
import pytest
from codedcal.workflow_validation import (ImageCase,Catalog,WorkflowResult,canonical,warp,apply_corrections,compare,assess,challenge_plan,legacy_challenge_plan,run_challenge)
from codedcal.workflow_bardensr import bounded_translation


def case(rounds=4,colors=4):
    rng=np.random.default_rng(11);f=rounds*colors
    im=rng.uniform(size=(f,1,40,40)).astype(np.float32)
    b=np.zeros((f,8))
    for k in range(8):b[np.arange(rounds)*colors+rng.integers(0,colors,size=rounds),k]=1
    return ImageCase(im,np.ones_like(im),b,np.repeat(np.arange(rounds),colors),np.tile(np.arange(colors),rounds),np.arange(8),np.array([7]))


def test_comparator_detects_count_preserving_label_swap_with_zero_composition_error():
    a=Catalog([[0,0,0],[0,0,10]],[0,1]);b=Catalog([[0,0,0],[0,0,10]],[1,0])
    d=compare(a,b,radius=1);assert d['identity_loss']==1 and d['composition_tv']==0
    assert not compare(Catalog([],[]),Catalog([],[]))['supported']


def test_empty_prediction_and_unsupported_reference_are_different():
    a=Catalog([[0,1,1]],[0]);r=WorkflowResult(Catalog([],[]),np.zeros((2,3)),1)
    assert assess(a,r,[],minimum_reference=1)['failed']
    assert assess(Catalog([],[]),r,[])['status']=='unsupported'


def test_native_inverse_and_gauge_have_exact_common_context():
    c=case();t=np.zeros((16,3));t[:,1:]=np.tile([[2,1],[-2,-1]],(8,1))
    obs=c.replace_images(warp(c.images,t),warp(c.support,t));restored,m=apply_corrections(obs,-t)
    expected=warp(c.images,np.repeat(t[:1],16,axis=0))
    np.testing.assert_array_equal(restored.images[:,:,8:-8,8:-8],expected[:,:,8:-8,8:-8])
    np.testing.assert_array_equal(canonical(-t+3),canonical(-t))
    assert np.all(m[0]==0)


def test_fractional_oracle_is_composition_not_false_exact_inverse():
    c=case();t=np.zeros((16,3));t[:,1]=np.tile([.25,-.25],8)
    obs=c.replace_images(warp(c.images,t),warp(c.support,t));restored,_=apply_corrections(obs,-t)
    assert not np.array_equal(restored.images[:,:,8:-8,8:-8],warp(c.images,np.repeat(t[:1],16,axis=0))[:,:,8:-8,8:-8])
    np.testing.assert_array_equal(restored.images,warp(warp(c.images,t),canonical(-t)))


@pytest.mark.parametrize('r,c',[(4,4),(10,4),(8,2),(6,3)])
def test_legacy_design_budget_holdout_and_common_histogram(r,c):
    a=case(r,c);p=legacy_challenge_plan(a)
    assert all(len(v)==4 for v in p['strategies'].values()) and len(p['heldout'])==5
    assert len(p['fractional'])==2 and p['logical_executions_per_workflow']==47
    assert not set(p['heldout'])&set(sum(p['strategies'].values(),[]))
    bank=p['selection']['bank'];signatures={tuple(sorted(x['mask'])) for x in bank};assert len(signatures)==1
    assert p==legacy_challenge_plan(a)


def test_bounded_correlation_reports_correction_sign():
    a=np.zeros((1,40,40));a[0,18,17]=1;b=warp(a[None],np.array([[0,2,-1]]))[0]
    result=bounded_translation(a,b,(0,5,5));np.testing.assert_allclose(result,[0,-2,1],atol=1e-8)


def test_adapter_does_not_receive_truth_or_reference_and_both_calls_charged():
    c=case();points=Catalog([[0,12+i,12] for i in range(12)],np.zeros(12,int))
    class Adapter:
        def __init__(self):self.n=0
        def run(self,image_case,policy):
            assert isinstance(image_case,ImageCase) and not hasattr(image_case,'imposed')
            self.n+=1;return WorkflowResult(points,np.zeros((16,3)),1)
    a=Adapter();ch=challenge_plan(c)['challenges'][0];row=run_challenge(c,a,'x',ch)
    assert a.n==row['actual_workflow_executions']==row['logical_workflow_executions']==2
    assert row['assessment']['status']=='scored'


def test_padded_integer_inverse_is_exact_globally_and_context_is_preserved():
    from codedcal.workflow_validation import pad_for_challenges
    c=case();p=pad_for_challenges(c)
    assert p.metadata['calibration_region']==[[0,1],[16,56],[16,56]]
    for ch in challenge_plan(c)['challenges']:
        if ch['role']=='fractional':continue
        t=np.array(ch['maps']);assert np.all(t[0]==0)
        restored=warp(warp(p.images,t),canonical(-t))
        np.testing.assert_array_equal(restored,p.images)
    assert pad_for_challenges(p) is p


def test_no_silent_fractional_class_or_readout_ids():
    c=case();c.rounds=c.rounds.astype(float);c.rounds[0]=.2
    with pytest.raises(ValueError,match='integers'):c.__post_init__()
    with pytest.raises(ValueError,match='integers'):Catalog([[0,1,2]],[.1])


def test_reject_all_has_no_adequate_acceptance_and_missing_is_not_success():
    from codedcal.workflow_decisions import bank_decision,confusion
    assert bank_decision([{'assessment':{'failed':None}}])['decision']=='unsupported'
    assert bank_decision([{'assessment':{'failed':True}},{'assessment':{'failed':None}}])['decision']=='reject'
    rows=[dict(truth={'decision':t},prediction={'decision':'reject'}) for t in ['accept','reject']]
    c=confusion(rows,'prediction');assert c['false_alert_rate']==1 and c['adequate_acceptance_rate']==0


def test_training_does_not_use_test_preparation_and_requires_adequate_controls():
    from codedcal.workflow_decisions import fit_threshold,leave_preparation_out
    rows=[dict(preparation='A',truth={'decision':'accept'},score=float(i)) for i in range(5)]
    rows += [dict(preparation='B',truth={'decision':'reject'},score=20.)]
    p=leave_preparation_out(rows,'score');assert not p[0]['rule']['supported']
    assert p[-1]['rule']['supported'] and p[-1]['rule']['adequate_training']==5 and p[-1]['decision']=='reject'
    assert not fit_threshold(rows[:4],'score')['supported']


def test_combined_qc_does_not_convert_missing_blanks_to_safety():
    from codedcal.workflow_decisions import fit_qc_union,qc_union_lopo
    rows=[dict(preparation='A',truth={'decision':'accept'},native_blank_proxy=.01,native_fit_loss=.1,native_inverse_yield=.01) for _ in range(5)]
    rows.append(dict(preparation='A',truth={'decision':'reject'},native_blank_proxy=.5,native_fit_loss=.1,native_inverse_yield=.01))
    rule=fit_qc_union(rows);assert rule['supported'] and rule['training_failure_detection']==1
    rows.append(dict(preparation='B',truth={'decision':'accept'},native_blank_proxy=None,native_fit_loss=.1,native_inverse_yield=.01))
    assert qc_union_lopo(rows)[-1]['decision']=='unsupported'


def test_provenance_reference_paths_are_removed_at_adapter_boundary():
    c=case();c.metadata['historical_reference']='/do/not/send/reference.json'
    class Adapter:
        def run(self,im,policy):
            assert 'historical_reference' not in im.metadata
            return WorkflowResult(Catalog([[0,30,30]]*12,[0]*12),np.zeros((len(im.images),3)),1)
    run_challenge(c,Adapter(),'x',challenge_plan(c)['challenges'][0])


def test_cli_plan_only_does_not_load_optional_adapter_or_overwrite(tmp_path):
    from codedcal.workflow_cli import main
    c=case();source=tmp_path/'input.npz';dest=tmp_path/'out'
    np.savez(source,images=c.images,support=c.support,codebook=c.codebook,rounds=c.rounds,colors=c.colors,class_ids=c.class_ids,blank_classes=c.blank_classes)
    assert main([str(source),'--output',str(dest),'--adapter','absent_module:absent_factory'])==0
    assert (dest/'plan.json').exists()
    with pytest.raises(FileExistsError):main([str(source),'--output',str(dest)])
    with pytest.raises(SystemExit):main([str(source),'--output',str(tmp_path/'other'),'--max-seconds','nan'])
