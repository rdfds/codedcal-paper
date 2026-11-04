"""Scientific controls for the new sensitivity, without touching old outputs."""
import numpy as np
from sensitivity import (settings, registry, event_phase_sums, phase_sums, O,
                         fast_clusters, fieller_set, signed_t, Q)
from response_uncertainty import rectangle_image, fit_response


def test_original_settings_are_historical():
    for row in registry():
        grid = settings(row)
        assert len(grid) == 9 and len({a['name'] for a in grid}) == 9
        base = grid[4]
        assert base['name'] == 'r1-original'
        assert base['radius_nm'] == (8. if row['tag'].startswith('sdc_') else 10.)
        assert base['qc']['min_per_panel'] == (3 if row['tag'].startswith('sdc_') else 2)
        assert base['qc']['max_zsd'] == (22. if row['tag'] == 'npc' else None)


def test_event_localization_membership_keeps_boundary_rows():
    # Mean event and individual row can occupy different panels or a guard gap.
    xy = np.array([[0., 1.], [2., 3.], [4., 5.], [6., 7.], [9., 8.]])
    panel = np.array([0, -1, 1, 2, 3]); event = np.array([0, 0, 0, 1, 1])
    z, n = event_phase_sums(xy, panel, event, 2)
    for chosen in ([0], [1], [0, 1], []):
        m = np.isin(event, chosen)
        expected_z, expected_n = phase_sums(xy[m], panel[m])
        assert np.allclose(z[chosen].sum(0), expected_z, atol=1e-14)
        assert np.array_equal(n[chosen].sum(0), expected_n)


def test_partitions_at_all_grid_radii_and_strict_boundary():
    rng = np.random.default_rng(2610231)
    for radius in (6.4, 8., 9.6, 10., 12.):
        for pos in (np.empty((0, 3)), np.array([[0, 0, 0], [radius, 0, 0], [radius*1.9, 0, 0]]),
                    rng.normal(size=(80, 3))*15):
            a = O.cluster_events(pos, thresh=radius); b = fast_clusters(pos, thresh=radius)
            assert np.array_equal(a[:, None] == a[None, :], b[:, None] == b[None, :])


def test_signed_equivalent_targets_and_empty_support():
    comp = np.ones((2, 2, 4)); comp[:, :, 0] = np.exp(Q[:, None]**2*2.)
    assert np.allclose(signed_t(comp), -2.)
    assert np.isnan(signed_t(np.zeros_like(comp))).all()
    assert fieller_set(np.zeros(12), np.zeros(12))['kind'] == 'all_real'


def test_rectangle_preserves_unsupported_and_negative_targets():
    assert rectangle_image([1, 2], [-.1, .1], [-1, 2])['interval'] is None
    assert rectangle_image(None, [-.1, .1], [1, 2])['interval'] is None
    result = rectangle_image([-2, -1], [0, 1], [1, 2])
    assert result['interval'] == [-3., -.5]
    # Coefficient uncertainty cannot yield a narrower set than plugging in
    # coefficient values contained by the rectangle.
    result = rectangle_image([1, 2], [-.1, .1], [.5, 2])
    assert result['interval'][0] <= 1 and result['interval'][1] >= 2


def test_existing_response_coefficients_and_intervals_replay():
    import json
    from sensitivity import OLD, WEIGHTS, student
    comparison = json.loads((OLD/'comparison.json').read_text())
    count = 0
    for acquisition in comparison['acquisitions']:
        path = OLD/f"controls-{acquisition['tag']}.json"
        if not path.exists():
            continue
        controls = json.loads(path.read_text())
        for j, name in enumerate(WEIGHTS):
            for model in ('gaussian', 'residual'):
                rows = [r for r in controls['records'] if r['model'] == model]
                beta, cov, _ = fit_response(rows, j)
                old = acquisition['response_sensitivity'][name]['models'][model]
                assert np.allclose(beta, [old['intercept'], old['slope']], atol=1e-10, rtol=0)
                interval = beta[1]+np.array([-1, 1])*student.ppf(.975, 3999)*np.sqrt(cov[1, 1])
                assert np.allclose(interval, old['slope_interval'], atol=1e-10, rtol=0)
                assert np.linalg.eigvalsh(cov).min() > -1e-12
                count += 1
    assert count == 36
