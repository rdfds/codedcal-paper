"""Design failures must be visible before any expensive workflow call."""
import copy
import json
import numpy as np
import pytest
from codedcal import workflow_cli
from codedcal import workflow_validation as wf
from codedcal.workflow_design import audit_plan,require_design


def case(rounds=4,colors=4):
    rng=np.random.default_rng(17);f=rounds*colors
    images=rng.uniform(size=(f,1,40,40)).astype(np.float32)
    book=np.zeros((f,8))
    for i in range(8):book[np.arange(rounds)*colors+rng.integers(colors,size=rounds),i]=1
    return wf.ImageCase(images,np.ones_like(images),book,np.repeat(np.arange(rounds),colors),np.tile(np.arange(colors),rounds),np.arange(8),np.array([7]))


@pytest.mark.parametrize('rounds,colors',[(4,4),(10,4),(8,2),(4,2),(6,4)])
@pytest.mark.parametrize('strategy',['mixed_conventional','mechanism_weighted'])
def test_each_family_and_independent_nonadditive_probe_are_required(rounds,colors,strategy):
    c=case(rounds,colors);p=wf.challenge_plan(c,strategy=strategy)
    assert p['design_status']=='supported' and p['design_audit']['supported']
    assert [r['mechanism_family'] for r in p['challenges']]==['round','color','interaction','frame']
    assert len(p['challenges'])==4 and p['logical_executions_per_workflow']==9
    assert p['heldout']==p['fractional']==[]
    assert p['design_audit']['strategies'][strategy]['sampled_nonadditive_rank']==2
    assert p==wf.challenge_plan(c,strategy=strategy)
    assert not p['selection']['weighted_variant_comparatively_validated']


def test_pair_coverage_does_not_satisfy_mechanism_coverage():
    c=case();old=wf.legacy_challenge_plan(c);bank=old['selection']['bank']
    selection=old['selection']['strategies']['readout_aware'];ids=old['strategies']['readout_aware']
    for token,index in zip(ids,selection):
        next(r for r in old['challenges'] if r['id']==token)['mechanism_family']=bank[index]['family']
    old['strategies']={'readout_aware':ids}
    audit=audit_plan(c,old)
    assert not audit['supported'] and audit['strategies']['readout_aware']['sampled_nonadditive_rank']==0
    assert any('interaction' in reason for reason in audit['reasons'])


def test_fake_family_labels_cannot_bypass_actual_map_guard():
    c=case();p=wf.challenge_plan(c)
    p['challenges'][2]['maps']=copy.deepcopy(p['challenges'][0]['maps'])
    with pytest.raises(ValueError,match='invalid_mechanism_maps:interaction'):require_design(c,p)


def test_repeated_interaction_is_not_an_independent_frame_probe():
    c=case();p=wf.challenge_plan(c)
    p['challenges'][3]['maps']=copy.deepcopy(p['challenges'][2]['maps'])
    with pytest.raises(ValueError,match='invalid_mechanism_maps:frame'):require_design(c,p)


def test_unsupported_design_spends_no_adapter_calls(tmp_path,monkeypatch):
    c=case(2,2);source=tmp_path/'input.npz'
    np.savez(source,images=c.images,support=c.support,codebook=c.codebook,rounds=c.rounds,colors=c.colors,class_ids=c.class_ids,blank_classes=c.blank_classes)
    def forbidden(*args):raise AssertionError('Must stop before loading the adapter')
    monkeypatch.setattr(workflow_cli,'import_adapter',forbidden)
    assert workflow_cli.main([str(source),'--output',str(tmp_path/'out'),'--run'])==2
    audit=json.loads((tmp_path/'out/design.json').read_text())
    assert not audit['supported'] and not (tmp_path/'out/qualification.json').exists()


def test_optional_research_panels_are_explicit_and_disjoint():
    c=case();p=wf.challenge_plan(c,include_holdout=True,include_fractional=True)
    assert p['design_status']=='supported' and len(p['heldout'])==5 and len(p['fractional'])==2
    assert len(p['challenges'])==11 and p['logical_executions_per_workflow']==23
    assert not set(p['heldout'])&set(p['strategies']['mixed_conventional'])


def test_run_api_rejects_incomplete_design_before_native_execution(tmp_path):
    c=case();p=wf.challenge_plan(c);p['strategies']['mixed_conventional']=p['strategies']['mixed_conventional'][:2]
    class Never:
        def run(self,*args):raise AssertionError('An incomplete design must not reach the adapter')
    with pytest.raises(ValueError,match='Unsupported challenge design'):wf.run_suite(c,Never(),['x'],tmp_path/'out',plan=p)
    assert not (tmp_path/'out').exists()
