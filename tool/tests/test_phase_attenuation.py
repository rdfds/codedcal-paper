"""Deterministic identity/support controls; no new empirical or outcome study."""
import itertools
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import t

from codedcal.phase_attenuation import (
    _modulus_image, _root_product, phase_components, phase_ratio_set,
    phase_attenuation_uncertainty, phase_budget,
)
from codedcal.smlm import components
from codedcal.uncertainty import _interval

ROOT = Path(__file__).resolve().parents[2]


def population_statistics(m, c):
    # Unequal positive reference masses do not need to match between windows.
    ra, rb = m*np.conj(c), np.conj(m)/np.conj(c)
    return np.array([[ra.real*2, ra.imag*2, 2.],
                     [rb.real*7, rb.imag*7, 7.]])


def test_retained_five_nm_mean_shift_witness_is_repaired_without_changing_old_result():
    path = ROOT/'outputs/full-evidence-audit-sep26/geometric-scope.json'
    witness = json.loads(path.read_text())
    q = np.array(witness['q']); true_T = witness['true_T']
    group = np.stack([population_statistics(np.exp(-qq**2*true_T), np.exp(1j*qq*5))
                      for qq in q])
    report = phase_attenuation_uncertainty(np.broadcast_to(group, (32, *group.shape)), q)
    assert report['all_frequencies_usable'] and not report['physical_assumptions_verified']
    for row in report['per_frequency']:
        interval = row['gaussian_equivalent_T_intervals'][0]
        assert interval['lower'] == pytest.approx(true_T)
        assert interval['upper'] == pytest.approx(true_T)
    # Preserve the actual historical failure; do not relabel it a success.
    np.testing.assert_allclose(witness['cases']['complex_mean_shift']['geometric_T'],
                               [22.058424044372288, 22.884439906626888])


def test_real_projection_is_nonidentifying_at_one_frequency_but_phase_separates_models():
    # Two valid population models have the SAME real C/R columns. One has
    # variance 9 + a kernel mean shift; the other has variance 22.058... and
    # no kernel change. This is the retained q=.1 witness, not a new dataset.
    witness=json.loads((ROOT/'outputs/full-evidence-audit-sep26/geometric-scope.json').read_text())
    q=witness['q'][0]; original_T=witness['true_T']
    alternative_T=witness['cases']['complex_mean_shift']['geometric_T'][0]
    left=population_statistics(np.exp(-q*q*original_T),np.exp(1j*q*5))
    right=population_statistics(np.exp(-q*q*alternative_T),1.)
    np.testing.assert_allclose(left[:,[0,2]],right[:,[0,2]],atol=1e-14)
    result=[phase_ratio_set(np.broadcast_to(x,(16,2,3)))['ratio'] for x in [left,right]]
    np.testing.assert_allclose(-np.log(result)/(q*q),[original_T,alternative_T],atol=1e-12)


@pytest.mark.parametrize('phase', [0., .5, np.pi/2, np.pi, -2.4])
@pytest.mark.parametrize('m', [.7, .7*np.exp(.8j)])
def test_population_identity_handles_complex_kernel_and_nonzero_motion_mean(phase, m):
    row = phase_ratio_set(np.broadcast_to(population_statistics(m, .8*np.exp(1j*phase)), (16,2,3)))
    assert row['usable'] and row['ratio'] == pytest.approx(abs(m))
    np.testing.assert_allclose([row['intervals'][0][key] for key in ('lower','upper')], abs(m))


def test_exact_code_and_motion_enumeration_preserves_weights_and_phase():
    # Three sources: enumerate all 4^3 permanent codes and 2^3 displacements.
    # This is a finite expectation identity, not 512 independent realizations.
    counts = np.array([[2,3,1,4],[4,1,2,3],[1,5,3,2]])
    theta = np.array([.3, 1.1, 2.4]); psi = np.array([.7,.9,.5])*np.exp(1j*np.array([.2,-.4,.6]))
    factor = .8*np.exp(.9j); angle = np.arccos(.75); means = []
    for signs in itertools.product((-1,1), repeat=3):
        displacement = np.array(signs)*angle+.4
        source = counts*psi[:,None]*np.exp(1j*theta[:,None])
        source[:,2:] *= factor*np.exp(1j*displacement[:,None])
        z = np.zeros((64,4,4),complex); n = np.zeros((64,4,4))
        for g, codes in enumerate(itertools.product(range(4),repeat=3)):
            for s, code in enumerate(codes):
                z[g,code] += source[s]; n[g,code] += counts[s]
        means.append(phase_components(z,n).mean(axis=0))
    average = np.mean(means,axis=0)
    row = phase_ratio_set(np.broadcast_to(average,(16,2,3)))
    assert row['ratio'] == pytest.approx(.75,abs=1e-13)
    assert average[0,2] != average[1,2]  # Do not assume equal window weights.


def test_real_component_agrees_with_unchanged_historical_aggregation():
    rng = np.random.default_rng(24)
    z = rng.normal(size=(16,4,4))+1j*rng.normal(size=(16,4,4))
    n = rng.uniform(.5,4,(16,4,4))
    old = components(z,n); new = phase_components(z,n)
    np.testing.assert_allclose(new[...,0],old[...,0],rtol=1e-13,atol=1e-13)
    np.testing.assert_allclose(new[...,2],old[...,1],rtol=1e-13,atol=1e-13)
    assert np.max(np.abs(new[...,1])) > 0


def test_fully_dependent_components_match_student_rectangle():
    x = np.array([.6,.7,.8,.9,1.])
    group = np.stack([np.column_stack([x*0,x,np.ones(5)]),
                      np.column_stack([x*0,-x,np.ones(5)])],axis=1)
    row = phase_ratio_set(group)
    critical = t.ppf((1+.9875)/2,4); radius = critical*x.std(ddof=1)/np.sqrt(5)
    assert row['component_confidence_level'] == .9875 and row['usable']
    np.testing.assert_allclose([row['intervals'][0][key] for key in ('lower','upper')],
                               [x.mean()-radius,x.mean()+radius])
    assert row['windows'][0]['real']['ratio'] == 0
    assert row['windows'][1]['imag']['ratio'] < 0  # Negative components are valid.
    assert row['windows'][0]['real']['usable'] and row['windows'][1]['imag']['usable']


def test_interval_propagation_keeps_weak_denominators_and_zero_support():
    reference = np.array([-.9,-.9,1.,1.])
    group = np.stack([np.column_stack([np.ones(4),np.ones(4),reference])]*2,axis=1)
    row = phase_ratio_set(group)
    assert all(w['real']['kind']=='disjoint' for w in row['windows'])
    assert row['intervals'][0]['upper'] is None and not row['usable']
    assert row['reason'] == 'unbounded_component_confidence_set'
    empty = phase_ratio_set(np.zeros((8,2,3)))
    assert empty['ratio'] is None and not empty['usable']
    assert empty['kind']=='unbounded_nonnegative'
    json.dumps(row,allow_nan=False);json.dumps(empty,allow_nan=False)


def test_disjoint_complex_rectangles_do_not_become_finite_by_clipping():
    parts = [_interval(upper=-2.,upper_closed=True), _interval(3.,4.,True,True)]
    image = _modulus_image(parts,[_interval(0.,0.,True,True)])
    assert image == [_interval(2.,None,True,False)]
    zero = [_interval(0.,0.,True,True)]
    assert _root_product(zero,[_interval(0.,None,False,False)]) == zero
    assert _modulus_image(None,zero) is None


def test_above_one_and_zero_magnitude_keep_their_actual_support():
    group = np.broadcast_to(population_statistics(1.2,.8j),(16,2,3))
    row = phase_attenuation_uncertainty(group[:,None],[.1])['per_frequency'][0]
    assert row['usable'] and row['unit_attenuation_excluded']
    assert not row['attenuation_model_compatible']
    assert row['gaussian_equivalent_T_intervals'][0]['upper'] < 0
    zero = phase_ratio_set(np.broadcast_to([[0.,0.,1.],[0.,0.,1.]],(8,2,3)))
    assert zero['ratio']==0 and not zero['usable']
    assert zero['reason']=='magnitude_not_separated_from_zero'


def test_source_dependent_kernel_change_remains_a_counterexample():
    a = .8; kernel = np.array([.8,.4])**2; change = np.array([.5,1.2])
    wa, wb = np.array([.9,.1]),np.array([.1,.9])
    ca,ra = a*np.sum(wa*kernel*change),np.sum(wa*kernel)
    cb,rb = a*np.sum(wb*kernel*change),np.sum(wb*kernel*change**2)
    row = phase_ratio_set(np.broadcast_to([[ca,0.,ra],[cb,0.,rb]],(16,2,3)))
    assert abs(row['ratio']-a)>.05  # A finite set is not a physical-model check.


@pytest.mark.parametrize('scale',[1e-180,1e180])
def test_common_scaling_preserves_complex_confidence_support(scale):
    x = np.linspace(.7,.9,32)
    group = np.stack([np.column_stack([x,.2*x,np.ones(32)]),
                      np.column_stack([x,-.2*x,np.ones(32)])],axis=1)
    before,after = phase_ratio_set(group),phase_ratio_set(group*scale)
    assert before['kind']==after['kind'] and before['usable']==after['usable']
    np.testing.assert_allclose([before['intervals'][0][key] for key in ('lower','upper')],
                               [after['intervals'][0][key] for key in ('lower','upper')])


@pytest.mark.parametrize('data',[np.ones((1,2,3)),np.full((8,2,3),np.nan)])
def test_invalid_component_sets_remain_explicit_and_json_safe(data):
    row = phase_ratio_set(data)
    assert row['kind']=='invalid' and row['intervals'] is None and not row['usable']
    json.dumps(row,allow_nan=False)


@pytest.mark.parametrize('level',[0.,1.,np.nan,np.nextafter(1.,0.)])
def test_invalid_confidence_rejected(level):
    with pytest.raises(ValueError,match='confidence_level'):
        phase_ratio_set(np.ones((8,2,3)),level)


def test_shape_frequency_and_count_validation():
    with pytest.raises(ValueError,match='shape'):
        phase_ratio_set(np.ones((8,2,2)))
    with pytest.raises(ValueError,match='real columns'):
        phase_ratio_set(np.ones((8,2,3),complex))
    for q in [[],[0.],[np.nan],[[.1]]]:
        with pytest.raises(ValueError,match='q must'):
            phase_budget(np.ones((8,4,4,1),complex),np.ones((8,4,4)),q)
    with pytest.raises(ValueError,match='nonnegative'):
        phase_components(np.ones((8,4,4),complex),-np.ones((8,4,4)))


@pytest.mark.parametrize('q',[1e-200,1e200])
def test_nonrepresentable_variance_transform_does_not_invent_a_statistical_infinity(q):
    group = np.broadcast_to(population_statistics(.8,.8j),(16,1,2,3))
    row = phase_attenuation_uncertainty(group,[q])['per_frequency'][0]
    assert row['usable'] and row['gaussian_equivalent_T_intervals'] is None
    assert row['gaussian_transform_reason']=='nonrepresentable_gaussian_transform'
    json.dumps(row,allow_nan=False)


@pytest.mark.parametrize('supported',[True,False])
def test_installed_command_path_preserves_phase_and_unresolved_support(tmp_path,capsys,supported):
    from codedcal import cli
    z = np.zeros((8,4,4,2),complex); n = np.zeros((8,4,4))
    q = np.array([.1,.15])
    if supported:
        z[:,0] = 1.; n[:,0] = 1.
        z[:,0,2:,:] = np.exp(1j*q*5)
    path = tmp_path/'input.npz';out = tmp_path/'report.json'
    np.savez(path,z=z,n=n,q=q)
    assert cli.main(['smlm-phase-budget',str(path),'--json',str(out)]) == (0 if supported else 2)
    result = json.loads(out.read_text(),parse_constant=lambda value:pytest.fail(value))
    assert result['phase_uncertainty']['all_frequencies_usable'] is supported
    for row in result['phase_uncertainty']['per_frequency']:
        if supported:
            assert row['ratio']==pytest.approx(1.)
            assert row['gaussian_equivalent_T_intervals'][0]['lower']==pytest.approx(0.,abs=1e-10)
        else:
            assert row['kind']=='unbounded_nonnegative' and row['intervals'][0]['upper'] is None
    text = capsys.readouterr().out
    assert 'Mean displacement and kernel phase are confounded' in text
    assert 'Physical assumptions are NOT verified' in text
    assert ('UNRESOLVED' in text) is not supported


def test_phase_command_reports_invalid_confidence_without_writing_a_result(tmp_path,capsys):
    from codedcal import cli
    path=tmp_path/'input.npz';np.savez(path,z=np.ones((8,4,4,1),complex),n=np.ones((8,4,4)),q=[.1])
    out=tmp_path/'must-not-exist.json'
    with pytest.raises(SystemExit) as error:
        cli.main(['smlm-phase-budget',str(path),'--confidence-level','1','--json',str(out)])
    assert error.value.code==2
    assert not out.exists() and '--confidence-level' in capsys.readouterr().err
