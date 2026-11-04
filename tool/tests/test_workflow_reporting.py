"""Scientific eligibility and end-to-end command controls, without native fits."""
import csv
import json
import numpy as np
import pytest

from codedcal import workflow_validation as wf
from codedcal import workflow_cli


def image_case():
    image=np.zeros((8,1,40,40),np.float32);image[:,:,20,20]=1
    book=np.array([[1,0],[0,1]]*4)
    return wf.ImageCase(image,np.ones_like(image),book,np.repeat(np.arange(4),2),np.tile([0,1],4),np.array([0,1]),np.array([],int))


def catalog(n=60):
    return wf.Catalog([[0,25+i/10,25] for i in range(n)],[0]*n)


class Adapter:
    def __init__(self):self.policies=[]
    def run(self,case,policy):
        self.policies.append(policy)
        return wf.WorkflowResult(catalog(0 if policy=='empty' else 60),np.zeros((len(case.images),3)),1,
            dict(actual_workflow_executions=1,actual_calibration_executions=1,actual_decoder_executions=1,nonblank_calls=0 if policy=='empty' else 60))


def test_reference_count_precedes_support_or_candidate_exception():
    c=image_case();small=catalog(4);poor=wf.WorkflowResult(catalog(),np.zeros((8,3)),.1)
    assert wf.assess(small,poor,[])['status']=='unsupported'
    class Broken:
        def __init__(self,n):self.calls=0;self.n=n
        def run(self,case,policy):
            self.calls+=1
            if self.calls==1:raise RuntimeError('candidate unavailable')
            return wf.WorkflowResult(catalog(self.n),np.zeros((8,3)),1)
    ch=wf.challenge_plan(c)['challenges'][0]
    unsupported=wf.run_challenge(c,Broken(4),'fixed',ch)['assessment']
    assert unsupported['failed'] is None and unsupported['reason']=='insufficient_reference_calls'
    failure=wf.run_challenge(c,Broken(60),'fixed',ch)['assessment']
    assert failure['failed'] is True and failure['reason']=='challenged_execution_failed'


def test_spatial_first_limit_does_not_invalidate_resolved_label_aware_match(monkeypatch):
    original=wf.spatial_match
    def bounded(a,b,radius):
        if len(a)>2:return np.empty((0,2),int),dict(unresolved_components=1)
        return original(a,b,radius)
    monkeypatch.setattr(wf,'spatial_match',bounded)
    a=wf.Catalog([[0,0,0],[0,0,1],[0,0,2],[0,0,3]],[0,0,1,1])
    result=wf.compare(a,a)
    assert result['supported'] and result['identity_loss']==0
    assert result['spatial_first_identity_loss'] is None and not result['spatial_first_supported']


def test_full_policy_native_controls_and_empty_native_gate():
    a=Adapter();gate=wf.qualify_input(image_case(),a,native_policy='full')
    assert gate['passed'] and a.policies==['full']*4 and gate['native_policy']=='full'
    assert not wf.qualify_input(image_case(),Adapter(),native_policy='empty')['passed']


def save_input(path):
    c=image_case()
    np.savez(path,images=c.images,support=c.support,codebook=c.codebook,rounds=c.rounds,colors=c.colors,class_ids=c.class_ids,blank_classes=c.blank_classes)


def test_cli_reports_controls_and_skips_unqualified_policy(tmp_path,monkeypatch):
    source=tmp_path/'input.npz';save_input(source);out=tmp_path/'result';adapter=Adapter()
    monkeypatch.setattr(workflow_cli,'import_adapter',lambda *args:adapter)
    assert workflow_cli.main([str(source),'--output',str(out),'--run','--native-policy','full','--policies','full,empty'])==0
    summary=json.loads((out/'summary.json').read_text());plan=json.loads((out/'plan.json').read_text())
    assert summary['stage']=='completed' and not summary['policies']['empty']['native_workflow_qualified']
    assert summary['policies']['empty']['heldout']['decision']=='unsupported'
    assert summary['logical_workflow_dispatches']==4+2+2*len(plan['challenges'])==len(adapter.policies)
    rows=list(csv.DictReader((out/'recovery.csv').open()))
    assert len([r for r in rows if r['role'].startswith('control:')])==4
    assert any(r['role']=='probe' and r['identity_loss']=='0.0' for r in rows)
    assert summary['policies']['full']['specified_challenges']['decision']=='all_specified_passed'
    assert summary['selected_strategy']=='mixed_conventional' and plan['heldout']==[]
    assert not list((out/'trials/empty').glob('*.json'))[1:]
    assert 'unknown-error risk' in (out/'summary.txt').read_text()


def test_cli_failed_controls_have_readable_evidence_and_no_trials(tmp_path,monkeypatch):
    source=tmp_path/'input.npz';save_input(source);out=tmp_path/'result'
    monkeypatch.setattr(workflow_cli,'import_adapter',lambda *args:Adapter())
    assert workflow_cli.main([str(source),'--output',str(out),'--run','--native-policy','empty'])==2
    assert json.loads((out/'summary.json').read_text())['stage']=='controls_failed'
    assert (out/'recovery.csv').exists() and not (out/'trials').exists()


def test_policy_gate_failure_is_not_a_qualified_outcome(tmp_path,monkeypatch):
    source=tmp_path/'input.npz';save_input(source);out=tmp_path/'result'
    monkeypatch.setattr(workflow_cli,'import_adapter',lambda *args:Adapter())
    assert workflow_cli.main([str(source),'--output',str(out),'--run','--native-policy','full','--policies','empty'])==2
    assert json.loads((out/'summary.json').read_text())['policies']['empty']['heldout']['decision']=='unsupported'


def test_deadline_retains_native_checkpoint_without_restarting(tmp_path,monkeypatch):
    now=[0.];monkeypatch.setattr(wf.time,'monotonic',lambda:now[0])
    class Slow(Adapter):
        def run(self,*args):now[0]+=2;return super().run(*args)
    with pytest.raises(TimeoutError):wf.run_suite(image_case(),Slow(),['full'],tmp_path/'trials',deadline_seconds=1)
    state=json.loads((tmp_path/'trials/state.json').read_text())
    assert state['stage']=='interrupted' and state['logical_workflow_dispatches']==1
    assert (tmp_path/'trials/full/native.json').exists()
    with pytest.raises(FileExistsError):wf.run_suite(image_case(),Adapter(),['full'],tmp_path/'trials')


def test_integer_oracle_mismatch_is_unsupported(tmp_path):
    class Nondeterministic(Adapter):
        def run(self,*args):
            r=super().run(*args)
            if len(self.policies)>1:r.catalog=wf.Catalog(r.catalog.points,[1]*60)
            return r
    c=image_case();plan=wf.challenge_plan(c)
    rows=wf.run_suite(c,Nondeterministic(),['full'],tmp_path/'trials',plan=plan)
    assert rows[0]['assessment']==dict(status='unsupported',reason='integer_oracle_native_mismatch',failed=None)
