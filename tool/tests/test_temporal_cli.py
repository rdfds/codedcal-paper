"""Exercise the real command, report and exit status on supported/weak data."""
import json

import numpy as np
import pytest

from codedcal import cli


def _input(tmp_path, supported):
    z = np.zeros((8, 4, 4, 1), complex)
    n = np.zeros((8, 4, 4))
    if supported:
        z[:, 0] = 1
        n[:, 0] = 1
    path = tmp_path/'events.npz'
    np.savez(path, z=z, n=n, q=np.array([.1]))
    return path


@pytest.mark.parametrize('supported', [True, False])
def test_command_status_and_strict_json(tmp_path, capsys, supported):
    path = _input(tmp_path, supported)
    report = tmp_path/'report.json'
    status = cli.main(['smlm-budget', str(path), '--json', str(report)])
    assert status == (0 if supported else 2)
    output = capsys.readouterr().out
    assert 'Fieller' in output and 'Physical assumptions are NOT verified' in output
    assert 'Legacy bounded fits' not in output
    data = json.loads(report.read_text(), parse_constant=lambda token: pytest.fail(token))
    row = data['uncertainty']['per_frequency'][0]
    assert row['kind'] == ('bounded' if supported else 'all_real')
    assert data['uncertainty']['all_frequencies_usable'] is supported
    assert 'T_pooled' in data and 'split' in data
    if not supported:
        assert 'UNRESOLVED' in output and '(-inf, +inf)' in output
        assert row['intervals'] == [dict(lower=None, upper=None,
                                          lower_closed=False, upper_closed=False)]


def test_legacy_fit_display_is_explicit_and_does_not_change_results(tmp_path, capsys):
    path = _input(tmp_path, True)
    report = tmp_path/'legacy.json'
    assert cli.main(['smlm-budget', str(path), '--legacy-fits',
                     '--confidence-level', '.9', '--json', str(report)]) == 0
    output = capsys.readouterr().out
    assert '90% Fieller' in output and 'Legacy bounded fits' in output
    assert json.loads(report.read_text())['uncertainty']['confidence_level'] == .9


@pytest.mark.parametrize('option,value', [('--q', '0'), ('--q', 'nan'),
    ('--q', 'inf'), ('--confidence-level', '1'), ('--confidence-level', 'nan')])
def test_invalid_options_fail_before_opening_input(option, value):
    with pytest.raises(SystemExit) as error:
        cli.main(['smlm-budget', 'does-not-exist.npz', option, value])
    assert error.value.code == 2


def test_inadequate_group_count_has_clear_error(tmp_path, capsys):
    path = tmp_path/'one.npz'
    np.savez(path, z=np.ones((1, 4, 4, 1), complex), n=np.ones((1, 4, 4)))
    assert cli.main(['smlm-budget', str(path), '--q', '.1']) == 2
    assert 'at least two independent groups' in capsys.readouterr().err


def test_negative_attenuation_is_retained_in_real_command(tmp_path, capsys):
    path = _input(tmp_path, True)
    with np.load(path) as data:
        z, n, q = data['z'], data['n'], data['q']
    z[:, 0, 2:] = -1  # within-window phase products positive; cross products negative
    np.savez(path, z=z, n=n, q=q)
    report = tmp_path/'negative.json'
    assert cli.main(['smlm-budget', str(path), '--json', str(report)]) == 2
    row = json.loads(report.read_text())['uncertainty']['per_frequency'][0]
    assert row['ratio'] == -1 and row['intervals'][0]['upper'] == -1
    assert row['gaussian_equivalent_T_intervals'] == []
    assert 'MODEL TENSION' in capsys.readouterr().out


def test_invalid_confidence_is_distinct_from_empty_set():
    assert cli._format_confidence_set(None) == 'unavailable'
    assert cli._format_confidence_set([]) == 'empty'


@pytest.mark.parametrize('variant', ['arithmetic', 'geometric', 'both'])
@pytest.mark.parametrize('supported', [True, False])
def test_requested_variants_define_status_and_both_are_retained_in_json(tmp_path, capsys, variant, supported):
    path = _input(tmp_path, supported)
    report = tmp_path/'variants.json'
    status = cli.main(['smlm-budget',str(path),'--uncertainty-variant',variant,
                       '--confidence-level','.9','--json',str(report)])
    assert status == (0 if supported else 2)
    data = json.loads(report.read_text(),parse_constant=lambda token: pytest.fail(token))
    expected = ['arithmetic','geometric'] if variant == 'both' else [variant]
    assert data['reported_uncertainty_variants'] == expected
    assert data['reported_uncertainty_usable'] is supported
    assert data['geometric_uncertainty']['window_confidence_level'] == .95
    assert data['uncertainty']['confidence_level'] == data['geometric_uncertainty']['confidence_level'] == .9
    output = capsys.readouterr().out
    assert ('geometric Fieller product sets' in output) == ('geometric' in expected)
    if 'geometric' in expected:
        assert 'window A:' in output and 'window B:' in output
        assert 'two 95% window sets' in output


def test_variants_can_disagree_without_automatic_switching(tmp_path, capsys):
    # Different window reference signs: pooled arithmetic is supported but
    # geometric support is forbidden. This exercises the actual aggregation.
    z = np.zeros((8,4,4,1),complex);n = np.zeros((8,4,4))
    n[:,0] = 1
    # Window references are 1 and -1/2, with cross moment 1/4:
    # pooled ratio = 1; window ratios = 1/4 and -1/2.
    z[:,0,:,0] = [1.,1.,1.,-.5]
    path = tmp_path/'different.npz';np.savez(path,z=z,n=n,q=[.1])
    default = tmp_path/'default.json';both = tmp_path/'both.json'
    assert cli.main(['smlm-budget',str(path),'--json',str(default)]) == 0
    assert cli.main(['smlm-budget',str(path),'--uncertainty-variant','both','--json',str(both)]) == 2
    a,b = json.loads(default.read_text()),json.loads(both.read_text())
    assert a['uncertainty']['all_frequencies_usable']
    assert not a['geometric_uncertainty']['all_frequencies_usable']
    assert a['uncertainty'] == b['uncertainty'] and a['geometric_uncertainty'] == b['geometric_uncertainty']
    assert 'No automatic variant selection' in capsys.readouterr().out


@pytest.mark.parametrize('column,value', [('group', .5), ('code', .5),
    ('panel', .5), ('x', float('nan')), ('weight', -1.)])
@pytest.mark.parametrize('format', ['csv', 'npz'])
def test_bad_events_cannot_silently_create_group_support(tmp_path, capsys, column, value, format):
    columns = ['group', 'code', 'panel', 'x', 'weight']
    events = np.array([[g, 0, p, 0, 1] for g in range(4) for p in range(4)], float)
    events[0, columns.index(column)] = value
    path = tmp_path/('events.'+format)
    if format == 'csv':
        np.savetxt(path, events, delimiter=',', header=','.join(columns), comments='')
    else:
        np.savez(path, events=events)
    assert cli.main(['smlm-budget', str(path)]) == 2
    assert 'smlm-budget:' in capsys.readouterr().err
