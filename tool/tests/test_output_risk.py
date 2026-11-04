"""Meaningful regressions for support, gauge and false-consensus risk."""
import importlib.util
from pathlib import Path

def test_output_risk_controls():
    root=Path(__file__).resolve().parents[2]
    spec=importlib.util.spec_from_file_location('risk_controls',root/'outputs/output-risk-oct06/check_controls.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    result=mod.run()
    assert result['passed'] and result['observable_false_consensus_flagged']
    assert result['invisible_label_permutation_flagged'] is False


def test_only_absent_codes_cannot_supply_decision_coverage():
    from codedcal import output_risk as O
    cat=dict(xyz=[[50,50,0]]*20,words=[99]*20,negative_ids=[99])
    comp=dict(abstained=False,compatible=[True])
    assert O.output_risk([cat],comp,0,[0])['abstained']
