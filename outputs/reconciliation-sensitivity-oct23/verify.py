"""Independent arithmetic/lineage recount of every sensitivity arm."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

P = Path(__file__).resolve().parent
ROOT = P.parents[1]
OLD = ROOT/'outputs/reference-reconciliation-oct07'
TAGS = ('npc', 'edf10rep1', 'edf10rep2', 'fig4untreated', 'rtx_rep1',
        'rtx_rep2', 'rtx_rep3', 'sdc_corrected', 'sdc_raw')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def equal(a, b):
    if a is None or b is None:
        return a is None and b is None
    return np.allclose(a, b, rtol=1e-10, atol=1e-8, equal_nan=True)


def main(partial=False, with_raw=False, write_receipt=False):
    checks = []; historical = json.loads((OLD/'comparison.json').read_text())
    narrowing_counts = {}
    dependencies = json.loads((P/'dependency-manifest.json').read_text())
    revisions = {r['path']: r for r in json.loads((P/'dependency-revisions.json').read_text())['revisions']}
    for dependency in dependencies['files']:
        revision = revisions.get(dependency['path'])
        if revision:
            assert dependency['sha256'] == revision['original_sha256'] == sha(ROOT/revision['original_copy'])
            assert sha(ROOT/dependency['path']) == revision['current_sha256']
        else:
            assert sha(ROOT/dependency['path']) == dependency['sha256'], dependency['path']
    additional = json.loads((P/'additional-dependencies.json').read_text())
    for dependency in additional['files']:
        assert sha(ROOT/dependency['path']) == dependency['sha256'], dependency['path']
    raw_inputs_checked = []
    tags = [tag for tag in TAGS if (P/f'result-{tag}.json').exists()] if partial else TAGS
    for tag in tags:
        r = json.loads((P/f'result-{tag}.json').read_text())
        d = np.load(P/f'groups-{tag}.npz'); old = np.load(OLD/f'reference-{tag}.npz')
        assert r['groups_sha256'] == sha(P/f'groups-{tag}.npz')
        assert r['old_reference_sha256'] == sha(OLD/f'reference-{tag}.npz')
        assert r['old_empirical_sha256'] == sha(OLD/f'empirical-{tag}.npz')
        assert equal(d['full_code'], old['full_code'])
        assert equal(d['selected_code'][4], old['selected_code'])
        assert equal(d['diagonal'][4], old['diagonal'])
        assert r['baseline_replay']['passed'] and r['baseline_replay']['source_rows'] == len(old['that'])
        assert len(r['arms']) == 9
        legacy = next(a for a in historical['acquisitions'] if a['tag'] == tag)
        # New sensitivity uses committed seed 26102361. Replaying old endpoint
        # values requires the old seed, not a different Monte Carlo sample.
        tile = d['tile']; ids = np.unique(tile)
        original_blocks = {key: np.stack([old[key][tile == k].sum(0) for k in ids])
            for key in ('selected_code', 'diagonal')}
        indices = np.random.default_rng(26100761).integers(len(ids), size=(2000, len(ids)))
        for new, old_pair in zip(r['arms'][4]['paired_selected_difference'],
                                 legacy['spatial_intervals']['paired_selected_difference']):
            assert equal(new['difference'], old_pair['code_minus_diagonal'])
            j = 0 if old_pair['q'] == .1 else 1
            ratios = []; denominators = []
            for key in ('selected_code', 'diagonal'):
                a = original_blocks[key][indices].sum(1)[:, j]
                denominator = a[:, :, 1].mean(1); denominators.append(denominator)
                ratios.append(a[:, :, 0].mean(1)/denominator)
            good = (ratios[0] > 0) & (ratios[1] > 0) & (denominators[0] > 0) & (denominators[1] > 0)
            difference = np.log(ratios[1][good]/ratios[0][good])/old_pair['q']**2
            assert equal(np.quantile(difference, [.025, .975]) if good.mean() >= .95 else None,
                old_pair['conditional_interval'])
        for m in r['input_members']:
            prior = next(v for v in legacy['input_members'] if v['round'] == m['round'])
            assert prior['sha256'] == m['sha256']
            if with_raw:
                # Resolve cache locations in this checkout; recorded absolute
                # paths are provenance, not portable input locations.
                import sys
                sys.path.insert(0, str(OLD))
                from reference_characterize import registry
                from common import input_path
                row = next(v for v in registry() if v['tag'] == tag)
                path = input_path(row, m['round'])
                assert sha(path) == m['sha256']
                raw_inputs_checked.append(str(path))
        for i, arm in enumerate(r['arms']):
            assert arm['sources'] == d['source_counts'][i]
            assert arm['localization_coverage'] == d['selected_localizations'][i]/d['all_localizations']
            assert arm['event_coverage'] == d['selected_events'][i]/d['all_events']
            assert equal(arm['original_row_jaccard'], d['overlap'][i, 0]/d['overlap'][i, 1]
                if d['overlap'][i, 1] else None)
            for key, array_key in [('selected_code', 'selected_code'), ('diagonal', 'diagonal')]:
                a = d[array_key][i].sum(0)
                ratio = (a[:, 0, 0]+a[:, 1, 0])/(a[:, 0, 1]+a[:, 1, 1])
                valid = (ratio > 0) & ((a[:, 0, 1]+a[:, 1, 1]) > 0)
                target = np.where(valid, -np.log(np.where(valid, ratio, 1.))/np.array([.1, .15])**2, np.nan)
                got = [v if v is not None else np.nan for v in arm[key]['point_T']]
                assert equal(target, got), (tag, i, key)
                # Every full confidence-set polynomial and topology is retained.
                assert len(arm[key]['uncertainty']['per_frequency']) == 2
            for j, name in enumerate(('unweighted', 'event_count', 'localization_count')):
                a = d['quadratic'][i, :, j].sum(0)
                assert equal(a[0]/a[1] if a[1] else None, arm['quadratic'][name]['point'])
                if i == 4:
                    assert equal(a[0]/a[1], legacy['estimates'][name]['raw'])
            for j, difference in enumerate(arm['paired_selected_difference']):
                left, right = arm['selected_code']['point_T'][j], arm['diagonal']['point_T'][j]
                assert equal(difference['difference'], left-right if left is not None and right is not None else None)
                assert 0 <= difference['valid_draws'] <= 2000
            expected = (arm['selected_to_diagonal_fit_gap'] < arm['full_to_diagonal_fit_gap']
                if arm['selected_to_diagonal_fit_gap'] is not None else None)
            assert arm['fit_gap_narrows'] == expected
        valid = [a for a in r['arms'] if all(description['fit'] is not None
            and not description['fit']['boundary']['arithmetic']
            and all(v is not None and np.isfinite(v) for v in description['point_T'])
            for description in (r['full_code'], a['selected_code'], a['diagonal']))]
        narrowing_counts[tag] = (sum(a['fit_gap_narrows'] for a in valid),
            sum(not a['fit_gap_narrows'] for a in valid), 9-len(valid))
        checks.append(dict(tag=tag, arms=9, baseline_replay=True, hashes=True, recounts=True))
    if partial:
        result = dict(stage='partial_recount_only', views=len(checks), arms=9*len(checks),
            response_propagation_verified=False, checks=checks)
        if write_receipt:
            (P/'verification-partial.json').write_text(json.dumps(result, indent=2)+'\n')
        print('Partial recount:', len(checks), 'views; no full-cohort conclusion')
        return
    summary = json.loads((P/'summary.json').read_text())
    assert summary['totals']['narrowed_arms'] == sum(v[0] for v in narrowing_counts.values())
    assert summary['totals']['worsened_arms'] == sum(v[1] for v in narrowing_counts.values())
    assert summary['totals']['unavailable_arms'] == sum(v[2] for v in narrowing_counts.values())
    for view in summary['views']:
        assert tuple(view[k] for k in ('narrowed', 'worsened', 'unavailable')) == narrowing_counts[view['tag']]
    # Separate diagnostic: identify deterministic single-source agreement and
    # independently fit the retained multi-source-only moment sums.
    from scipy.optimize import minimize_scalar
    secondary = json.loads((P/'singleton-diagnostic.json').read_text())
    secondary_counts = [0, 0, 0]
    for view in secondary['views']:
        d = np.load(P/f"groups-{view['tag']}.npz")
        for i, arm in enumerate(view['arms']):
            count = d['quadratic'][i, :, 0, 1]; single = count == 1; multi = count > 1
            assert arm['singleton_groups'] == int(single.sum())
            assert arm['multi_source_groups'] == int(multi.sum())
            assert equal(d['selected_code'][i, single], d['diagonal'][i, single])
            if not arm['supported']:
                secondary_counts[2] += 1; assert arm['gap_narrows'] is None; continue
            fits = []
            for array in (d['full_code'], d['selected_code'][i], d['diagonal'][i]):
                a = array[multi].sum(0)
                ratio = (a[:, 0, 0]+a[:, 1, 0])/(a[:, 0, 1]+a[:, 1, 1])
                fits.append(float(minimize_scalar(lambda t: np.sum((ratio-np.exp(-np.array([.1, .15])**2*t))**2),
                    bounds=(-30, 60), method='bounded').x))
            old_gap, new_gap = abs(fits[0]-fits[2]), abs(fits[1]-fits[2])
            assert abs(old_gap-arm['full_to_diagonal_fit_gap']) < 1e-7
            assert abs(new_gap-arm['selected_to_diagonal_fit_gap']) < 1e-7
            assert (new_gap < old_gap) == arm['gap_narrows']
            secondary_counts[0 if new_gap < old_gap else 1] += 1
    assert secondary_counts == [secondary['totals'][k] for k in ('narrowed', 'worsened', 'unavailable')]
    response = json.loads((P/'response-uncertainty.json').read_text())
    assert len(response['results']) == 9
    assert sum(r['status'] == 'conditional_propagation_completed' for r in response['results']) == 6
    images = 0
    for r in response['results']:
        if 'weights' not in r:
            assert r['tag'] in ('npc', 'sdc_corrected', 'sdc_raw'); continue
        for weight in r['weights'].values():
            f = weight['observation']
            obs = [f['intervals'][0]['lower'], f['intervals'][0]['upper']] if f['kind'] == 'bounded' else None
            for model in weight['models'].values():
                a, b = np.array(model['coefficient_intervals'])
                image = model['propagated']
                if b[0] <= 0 or obs is None:
                    assert image['interval'] is None
                else:
                    points = [(v-aa)/bb for v in obs for aa in a for bb in b]
                    assert equal([min(points), max(points)], image['interval'])
                images += 1
    artifacts = {p.name: sha(p) for p in sorted(P.iterdir()) if p.is_file()
        and p.name != 'verification.json' and not p.name.endswith('.log')}
    output = dict(passed=True, distinct_acquisitions=8, views=9, arms=81,
        response_images=images, dependency_hashes_checked=len(dependencies['files']),
        additional_dependency_hashes_checked=len(additional['files']),
        raw_input_checks=len(raw_inputs_checked), unique_raw_paths=len(set(raw_inputs_checked)),
        raw_inputs_freshly_verified=with_raw,
        singleton_diagnostic_recount=secondary_counts,
        checks=checks, artifact_sha256=artifacts,
        physical_validation=False, independent_confirmation=False,
        scope='Retained input identity and independent aggregation/rectangle recounts; not physical source truth.')
    if write_receipt:
        (P/'verification.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({k: v for k, v in output.items() if k not in ('checks', 'artifact_sha256')}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--partial', action='store_true')
    parser.add_argument('--with-raw', action='store_true', help='Also require original raw caches and freshly hash them')
    parser.add_argument('--write-receipt', action='store_true', help='Explicit maintainer mutation: refresh in-tree verification receipt; default is read-only')
    args = parser.parse_args(); main(args.partial, args.with_raw, args.write_receipt)
