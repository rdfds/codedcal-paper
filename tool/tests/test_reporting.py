"""Admission/support regressions; no new empirical observations or study."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from codedcal import smlm
from codedcal.reporting import budget


def aggregates():
    z = np.zeros((32, 4, 4, 1), complex)
    n = np.zeros((32, 4, 4))
    z[:, 0] = 1.; n[:, 0] = 1.
    z[:, 0, 2:, :] = np.exp(1j*np.pi/3)
    return z, n


def test_supported_report_preserves_every_historical_field():
    z, n = aggregates()
    with np.errstate(all='ignore'):
        original = smlm.budget(z, n, [.1])
    current = budget(z, n, [.1])
    guard = current.pop('numerical_reporting_guard')
    # Historical diagnostic NaNs compare unequal to themselves; their exact
    # serialized representations must still agree rather than being discarded.
    assert json.dumps(current, sort_keys=True) == json.dumps(original, sort_keys=True)
    assert guard['passed'] and not guard['historical_estimator_changed']
    assert not guard['physical_assumptions_verified']


@pytest.mark.parametrize('q', [1e-200, 1e-160, 1e160, 1e200])
def test_overflow_and_underflow_cannot_be_reported_as_bounds(q):
    with pytest.raises(ValueError, match='frequency squared.*rescale'):
        budget(*aggregates(), [q])


@pytest.mark.parametrize('q', [[.1+1j], [0.], [-1.], [np.inf], [np.nan]])
def test_invalid_or_complex_frequency_is_not_coerced(q):
    with pytest.raises(ValueError):
        budget(*aggregates(), q)


def test_coordinate_unit_rescaling_preserves_attenuation_and_scales_variance():
    # Aggregates depend on q*x, so changing x units and reciprocal q leaves
    # these phase/count inputs unchanged. T must change with squared units.
    z, n = aggregates()
    old = budget(z, n, [.1]); new = budget(z, n, [10.])
    for name in ('uncertainty', 'geometric_uncertainty'):
        left, right = [r[name]['per_frequency'][0] for r in (old, new)]
        assert left['intervals'] == right['intervals']
        for a, b in zip(left['gaussian_equivalent_T_intervals'], right['gaussian_equivalent_T_intervals']):
            for key in ('lower', 'upper'):
                assert b[key] == pytest.approx(a[key]/10000)


def test_complex_counts_rejected():
    z, n = aggregates()
    with pytest.raises(ValueError, match='counts must be real'):
        budget(z, n.astype(complex)+1j, [.1])


def test_zero_reference_remains_an_unresolved_statistical_result():
    z, n = aggregates()
    result = budget(z*0, n*0, [.1])
    assert not result['uncertainty']['all_frequencies_usable']
    assert not result['geometric_uncertainty']['all_frequencies_usable']
    assert result['numerical_reporting_guard']['passed']


def test_nonfinite_transform_is_rejected_even_inside_frequency_admission(monkeypatch):
    z, n = aggregates()
    with np.errstate(all='ignore'):
        report = smlm.budget(z, n, [.1])
    report['uncertainty']['per_frequency'][0]['gaussian_equivalent_T_intervals'][0]['lower'] = np.inf
    monkeypatch.setattr(smlm, 'budget', lambda *args: copy.deepcopy(report))
    with pytest.raises(ValueError, match='nonrepresentable Gaussian'):
        budget(z, n, [.1])


def test_npz_cli_rejects_extreme_units_without_writing_misleading_json(tmp_path):
    z, n = aggregates()
    path = tmp_path/'extreme.npz'; output = tmp_path/'report.json'
    np.savez(path, z=z, n=n, q=[1e160])
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]/'src'))
    run = subprocess.run([sys.executable, '-m', 'codedcal', 'smlm-budget', str(path),
                          '--json', str(output)], env=env, capture_output=True, text=True)
    assert run.returncode == 2 and 'rescale' in run.stderr
    assert not output.exists()


@pytest.mark.parametrize('command', ['smlm-budget', 'smlm-phase-budget'])
def test_npz_complex_frequency_is_not_discarded_by_either_command(tmp_path, command):
    z, n = aggregates()
    path = tmp_path/'complex.npz'; output = tmp_path/'report.json'
    np.savez(path, z=z, n=n, q=[.1+1j])
    env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parents[1]/'src'))
    run = subprocess.run([sys.executable, '-m', 'codedcal', command, str(path),
                          '--json', str(output)], env=env, capture_output=True, text=True)
    assert run.returncode == 2 and 'real' in run.stderr
    assert not output.exists()
