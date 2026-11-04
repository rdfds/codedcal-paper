"""Regression tests for user-visible search and QC failure contracts."""
import json
import numpy as np
import pytest
from codedcal import cli


def test_cli_search_recovers_motion_beyond_default_window(tmp_path, capsys):
    rng = np.random.default_rng(72)
    base = rng.uniform(100, 380, (800, 2))
    spots = np.vstack([np.column_stack([np.full(len(base), r), np.full(len(base), c),
                        base + [6.2*r, -4.1*r] + rng.normal(0, .05, base.shape)])
                        for r in range(8) for c in range(3)
                        if (r, c) not in ((6, 1), (7, 0))])
    path = tmp_path / 'field.npz'
    np.savez(path, spots=spots)
    report = tmp_path / 'trace.json'
    args = ['fish-trace', str(path), '--bins', '512', '--bin-width', '1',
            '--json', str(report)]
    assert cli.main(args) == 2
    assert 'SHIFT SEARCH TOO NARROW' in capsys.readouterr().out
    assert cli.main(args + ['--search', '80']) == 0
    result = json.loads(report.read_text())
    assert result['search_px'] == 80
    assert result['n_pairs_peak_outside_window'] == 0
    np.testing.assert_allclose(result['trace_signal_px'],
                               [[6.2*r, -4.1*r] for r in range(8)], atol=.15)


@pytest.mark.parametrize('radius', ['0', '-5', 'nan', 'inf'])
def test_invalid_search_radius_is_rejected(radius):
    with pytest.raises(SystemExit) as error:
        cli.main(['fish-trace', 'unused.npz', '--search', radius])
    assert error.value.code == 2


@pytest.mark.parametrize('passed', [True, False])
@pytest.mark.parametrize('round1', [True, False])
def test_magnification_qc_exit_status_and_report(tmp_path, monkeypatch, passed, round1):
    path = tmp_path / 'field.npz'
    np.savez(path, spots=np.array([[0, 0, 5, 5], [0, 1, 5, 5]]))
    estimate = dict(gate_pass=passed, gate_reason='test gate')
    monkeypatch.setattr(cli.chromatic, 'round1_qc', lambda *a, **kw: estimate)
    monkeypatch.setattr(cli.chromatic, 'estimate_scales', lambda *a, **kw: estimate)
    monkeypatch.setattr(cli.chromatic, 'gate', lambda *a: (passed, 'test gate'))
    report = tmp_path / 'qc.json'
    args = ['magnification', str(path), '--json', str(report)]
    assert cli.main(args + (['--round1'] if round1 else [])) == (0 if passed else 2)
    result = json.loads(report.read_text())
    assert result['gate_pass'] is passed
    assert 'not a confidence interval' in result['uncertainty_status']
    assert 'not an accuracy guarantee' in result['gate_interpretation']
    assert result['scale_se_gauge'] == '0'
