import hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent.parent.parent;sys.path.insert(0,str(ROOT/'tool/src'))
from codedcal import output_risk_weighted as W
HERE=Path(__file__).resolve().parent

def main():
    profile=dict(residual=np.zeros((3,16)),energy=np.ones((3,16)),block_counts=np.ones(16),coverage=1.)
    selection=W.selection_weights(profile);assert np.allclose(selection['weights'],1/3) and selection['multiple_tie_draws']==500
    base=dict(xyz=[[45+i,45,0] for i in range(16)],words=list(range(16)))
    changed=base|dict(words=list(range(1,16))+[0]);maps=np.zeros((3,8,3),int);risk=W.weighted_output_risk([base,changed,changed],selection,0,maps);assert abs(risk['identity_risk']-2/3)<1e-12 and risk['flagged'] and not risk['posterior_RNA_error_probability']
    superior=profile|dict(residual=np.r_[np.zeros((1,16)),np.ones((2,16))]);s=W.selection_weights(superior);assert s['weights']==[1.,0.,0.];safe=W.weighted_output_risk([base,changed,changed],s,0,maps);assert not safe['flagged']
    assert W.selection_weights(superior|dict(energy=superior['energy']*7,residual=superior['residual']*7))['weights']==s['weights']
    assert W.selection_weights(profile|dict(coverage=0.))['abstained']
    out=dict(passed=True,ambiguous_false_consensus_flags=True,ties_shared=True,dominated_candidates_do_not_force_warning=True,global_scale_invariant=True,missing_support_abstains=True,posterior_RNA_error_probability=False,script_sha256=hashlib.sha256(Path(W.__file__).read_bytes()).hexdigest(),protocol_sha256=hashlib.sha256((HERE/'WEIGHTED_REPAIR_PROTOCOL.md').read_bytes()).hexdigest());(HERE/'weighted-controls.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
