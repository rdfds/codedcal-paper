"""Reaggregate reporting details from retained reductions; no fitting or new draws.

Use --root for either the original archive or the independent research companion.
Output must be a new directory. Original qualification decisions stay unchanged.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

TAGS = ('npc', 'edf10rep1', 'edf10rep2', 'fig4untreated', 'rtx_rep1',
        'rtx_rep2', 'rtx_rep3', 'sdc_corrected', 'sdc_raw')


def export(root, out):
    out.mkdir(parents=True, exist_ok=False)
    pins = []

    def load(name, arrays=False):
        p = root/name
        pins.append(dict(path=name, bytes=p.stat().st_size,
                         sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        return np.load(p, allow_pickle=False) if arrays else json.loads(p.read_text())

    old = load('outputs/reference-reconciliation-oct07/comparison.json')['acquisitions']
    baselines, support, residuals = [], [], []
    for tag in TAGS:
        v = load(f'outputs/reconciliation-sensitivity-oct23/result-{tag}.json')
        with load(f'outputs/reconciliation-sensitivity-oct23/groups-{tag}.npz', True) as d:
            legacy = next(a for a in old if a['tag'] == tag)
            full = v['full_code']['fit']['total_arithmetic']
            base = v['arms'][4]
            assert base['name'] == 'r1-original'
            historical = legacy['historical_code_variants']['pooled']
            if tag == 'npc':
                historical = historical['all']
            assert full == legacy['full_code']['total_arithmetic']
            baselines.append(dict(view=tag,
                historical_full_phase_convention=historical['total_arithmetic'],
                full_centered_float64=full,
                selected_float64=base['selected_code']['fit']['total_arithmetic'],
                diagonal_float64=base['diagonal']['fit']['total_arithmetic']))
            tile_ids = np.unique(d['tile'])
            assert 1 <= len(tile_ids) <= 12 and set(tile_ids) <= set(range(12))
            for i, a in enumerate(v['arms']):
                counts = np.array([d['quadratic'][i, d['tile'] == t, 0, 1].sum()
                                   for t in tile_ids])
                assert counts.sum() == a['sources'] == d['source_counts'][i]
                assert np.all(counts == counts.astype(int))
                assert a['localization_coverage'] == d['selected_localizations'][i]/d['all_localizations']
                assert a['event_coverage'] == d['selected_events'][i]/d['all_events']
                gap = abs(full-a['diagonal']['fit']['total_arithmetic']) - abs(
                    a['selected_code']['fit']['total_arithmetic']-a['diagonal']['fit']['total_arithmetic'])
                assert np.isclose(gap, a['full_to_diagonal_fit_gap']-a['selected_to_diagonal_fit_gap'],
                                  rtol=1e-12, atol=1e-12)
                assert (gap > 0) == a['fit_gap_narrows']
                for j, p in enumerate(a['paired_selected_difference']):
                    support.append(dict(view=tag, setting=a['name'], q=p['q'],
                        sources=a['sources'], nonempty_blocks=int((counts > 0).sum()),
                        singleton_blocks=int((counts == 1).sum()), blocks=len(tile_ids),
                        partition_cells=12, resampling_tile_ids=tile_ids.tolist(),
                        source_counts_per_block=counts.astype(int).tolist(),
                        valid_draws=p['valid_draws'], total_draws=p['total_draws'],
                        valid_fraction=p['valid_draws']/p['total_draws'],
                        code_input_fraction=a['localization_coverage'],
                        event_fraction=a['event_coverage']))
            for j, q in enumerate((.1, .15)):
                kinds = dict(positive=0, negative=0, zero_spanning=0, unavailable=0)
                for a in v['arms']:
                    ci = a['paired_selected_difference'][j]['conditional_pointwise_interval']
                    kind = ('unavailable' if ci is None else 'positive' if ci[0] > 0
                            else 'negative' if ci[1] < 0 else 'zero_spanning')
                    kinds[kind] += 1
                assert sum(kinds.values()) == 9
                residuals.append(dict(view=tag, q=q, **kinds))

    rows = load('outputs/temporal-variants-oct17/summary.json')['summaries']
    gates = []
    for method in ('coded', 'known_source'):
        for variant in ('arithmetic', 'geometric'):
            chosen = [r for r in rows if r['primary'] and r['scope'] == 'observed'
                      and r['method'] == method and r['variant'] == variant]
            assert len(chosen) == 144
            failures = [set(r['failed_gates']) for r in chosen]
            assert all(r['qualified'] == (not f) for r, f in zip(chosen, failures))
            allowed = sum(not f or (f == {'nominal_coverage_in_bootstrap'} and
                r['coverage_bootstrap_nominal_failure_direction'] == 'above_nominal')
                for r, f in zip(chosen, failures))
            union = sum(bool(f & {'ratio_bias', 'usable_fraction', 'usable_coverage'}) for f in failures)
            gates.append(dict(method=method, variant=variant, conditions=144,
                all_inherited_gates=sum(r['qualified'] for r in chosen),
                descriptive_allow_above_nominal_only=allowed,
                bias_usability_or_joint_coverage_failure_union=union))

    minimum = min(support, key=lambda r: r['code_input_fraction'])
    summary = dict(views=9, sensitivity_settings=81, frequency_evaluations=len(support),
        minimum_nonempty_blocks=min(r['nonempty_blocks'] for r in support),
        minimum_valid_fraction=min(r['valid_fraction'] for r in support),
        minimum_coverage=minimum, qualification_thresholds_changed=False,
        new_fits=0, new_random_draws=0, source_pins=pins)
    for name, records in [('baselines', baselines), ('support', support),
                          ('residuals', residuals), ('variant-gates', gates)]:
        with (out/(name+'.csv')).open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=records[0], lineterminator='\n')
            w.writeheader(); w.writerows(records)
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    text = ['# Supplementary reporting details', '',
        'All entries below reaggregate retained reductions without new fits or random draws.', '',
        '## Full-population arithmetic baselines', '',
        'All fits have nm²-equivalent units. Historical phases used float32 for CD20 and SDC; NPC already used float64. The centered-float64 full-population baseline is fixed across the nine settings of a view and is exactly the full-code fit already plotted in Figure 2. Selected and diagonal columns use the original eligibility setting (r1-original); neither is physical ground truth. Thus Figure 2 already compares centered-float64 full and selected reductions; the separate earlier historical fits quantify the small arithmetic change, rather than a hidden confound within that display.', '',
        '| View | Full historical phase convention | Full centered float64 (Figure 2) | Selected float64 | Diagonal float64 |',
        '|---|---:|---:|---:|---:|']
    for r in baselines:
        text.append('| '+r['view']+' | '+' | '.join(f'{r[k]:.6f}' for k in (
            'historical_full_phase_convention','full_centered_float64','selected_float64','diagonal_float64'))+' |')
    text += ['', '## Spatial support and validity across every sensitivity setting', '',
        f"The minimum selected-source support is {summary['minimum_nonempty_blocks']} nonempty spatial blocks. The 4-by-3 partition has 12 possible cells; the actual resampling frame uses only tile IDs occupied by the original coarse groups (7–12 cells by view). Cells in that frame with zero selected sources remain in it. The minimum sign-valid draw fraction is {100*summary['minimum_valid_fraction']:.2f}% across all 162 frequency evaluations. A high valid fraction does not imply broad spatial support or calibrated coverage.", '',
        f"The minimum code-input coverage is {100*minimum['code_input_fraction']:.6f}% in {minimum['view']}, {minimum['setting']}, with {minimum['sources']} selected sources and {minimum['nonempty_blocks']} nonempty blocks. It rounds to 0.1% in the retained figure. Event coverage divides selected candidate binding-event rows by all candidate event rows in valid temporal panels; code-input coverage instead divides retained frame-level localization rows by all code-input localization rows in those panels. Burst multiplicity makes these denominators different.", '',
        'The complete support.csv contains block-level selected-source counts and valid-draw counts for all 81 settings at both frequencies. Three unsupported settings each evaluate two estimators at two frequencies: 3 × 2 × 2 = 12 all-real or disjoint individual Fieller sets.', '',
        '| View | q | Positive | Negative | Zero-spanning | Unavailable |',
        '|---|---:|---:|---:|---:|---:|']
    for r in residuals:
        text.append('| '+' | '.join(str(r[k]) for k in ('view','q','positive','negative','zero_spanning','unavailable'))+' |')
    text += ['', '## Temporal-variant qualification accounting', '',
        'The inherited four-gate qualification is unchanged. Its usable-and-covering frequency gate is at least 0.90, not 0.95. The final inherited gate asks whether the field-bootstrap range contains nominal 0.95. The companion column only ignores a failure of that last gate from above; it is a descriptive reclassification, not a new qualification criterion or a one-sided 95% confidence guarantee. The union counts overlapping bias, usability and joint-coverage failures once.', '',
        '| Method | Variant | All inherited gates /144 | Companion allowing above-nominal only /144 | Bias/usability/joint-coverage failure union /144 |',
        '|---|---|---:|---:|---:|']
    for r in gates:
        text.append('| '+' | '.join(str(r[k]) for k in ('method','variant','all_inherited_gates',
            'descriptive_allow_above_nominal_only','bias_usability_or_joint_coverage_failure_union'))+' |')
    (out/'report.md').write_text('\n'.join(text)+'\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'source_pins'}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); export(a.root.resolve(), a.output.resolve())
