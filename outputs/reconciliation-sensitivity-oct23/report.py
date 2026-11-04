"""Summarize every specified setting; never select a preferred agreement arm."""
import csv
from sensitivity import *


def fmt(value, digits=3):
    return 'unavailable' if value is None else f'{value:.{digits}f}'


def supported_fit_comparison(arm, full):
    descriptions = (full, arm['selected_code'], arm['diagonal'])
    return all(d['fit'] is not None and not d['fit']['boundary']['arithmetic']
        and all(v is not None and np.isfinite(v) for v in d['point_T']) for d in descriptions)


def run():
    views = [json.loads((P/f"result-{r['tag']}.json").read_text()) for r in registry()]
    responses = json.loads((P/'response-uncertainty.json').read_text())
    summary = []; csv_rows = []
    for view in views:
        arms = view['arms']; supported = [a for a in arms if supported_fit_comparison(a, view['full_code'])]
        ranges = {}
        for name in ('sources', 'localization_coverage', 'original_row_jaccard',
                     'full_to_diagonal_fit_gap', 'selected_to_diagonal_fit_gap'):
            values = [a[name] for a in arms if a[name] is not None]
            ranges[name] = [min(values), max(values)] if values else None
        residuals = []
        for j, q in enumerate(Q):
            positive = negative = unavailable = 0
            for a in arms:
                p = a['paired_selected_difference'][j]; ci = p['conditional_pointwise_interval']
                if ci is None:
                    unavailable += 1
                elif ci[0] > 0:
                    positive += 1
                elif ci[1] < 0:
                    negative += 1
            residuals.append(dict(q=float(q), positive_pointwise_intervals=positive,
                negative_pointwise_intervals=negative, unavailable_intervals=unavailable))
        result = dict(tag=view['tag'], narrowed=sum(a['fit_gap_narrows'] is True for a in supported),
            worsened=sum(a['fit_gap_narrows'] is False for a in supported),
            unavailable=9-len(supported), ranges=ranges, residuals=residuals,
            all_settings_narrow=bool(supported and len(supported) == 9 and
                all(a['fit_gap_narrows'] for a in supported)),
            all_settings_worsen=bool(supported and len(supported) == 9 and
                not any(a['fit_gap_narrows'] for a in supported)))
        result['arithmetic_boundary_arms'] = sum(any(
            a[key]['fit'] is not None and a[key]['fit']['boundary']['arithmetic']
            for key in ('selected_code', 'diagonal')) for a in arms)
        result['both_frequencies_fieller_usable_arms'] = sum(all(
            a[key]['uncertainty']['all_frequencies_usable']
            for key in ('selected_code', 'diagonal')) for a in arms)
        summary.append(result)
        for a in arms:
            row = dict(tag=view['tag'], setting=a['name'], sources=a['sources'],
                localization_coverage=a['localization_coverage'], event_coverage=a['event_coverage'],
                original_row_jaccard=a['original_row_jaccard'],
                full_to_diagonal_fit_gap=a['full_to_diagonal_fit_gap'],
                selected_to_diagonal_fit_gap=a['selected_to_diagonal_fit_gap'],
                fit_gap_narrows=a['fit_gap_narrows'])
            row['fit_comparison_supported'] = supported_fit_comparison(a, view['full_code'])
            row['selected_fit_boundary'] = (a['selected_code']['fit']['boundary']['arithmetic']
                if a['selected_code']['fit'] is not None else None)
            row['diagonal_fit_boundary'] = (a['diagonal']['fit']['boundary']['arithmetic']
                if a['diagonal']['fit'] is not None else None)
            for name in WEIGHTS:
                row[name+'_quadratic'] = a['quadratic'][name]['point']
            for j, q in enumerate(Q):
                paired_result = a['paired_selected_difference'][j]
                row[f'q{q:g}_selected_T'] = a['selected_code']['point_T'][j]
                row[f'q{q:g}_diagonal_T'] = a['diagonal']['point_T'][j]
                row[f'q{q:g}_difference'] = paired_result['difference']
                ci = paired_result['conditional_pointwise_interval']
                row[f'q{q:g}_lower'] = ci[0] if ci is not None else None
                row[f'q{q:g}_upper'] = ci[1] if ci is not None else None
                row[f'q{q:g}_valid_bootstraps'] = paired_result['valid_draws']
                row[f'q{q:g}_selected_set'] = a['selected_code']['uncertainty']['per_frequency'][j]['kind']
                row[f'q{q:g}_diagonal_set'] = a['diagonal']['uncertainty']['per_frequency'][j]['kind']
            csv_rows.append(row)
    with (P/'all-arms.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(csv_rows[0])); writer.writeheader(); writer.writerows(csv_rows)
    totals = dict(views=9, distinct_acquisitions=8, arms=81,
        all_settings_narrow=sum(s['all_settings_narrow'] for s in summary),
        all_settings_worsen=sum(s['all_settings_worsen'] for s in summary),
        narrowed_arms=sum(s['narrowed'] for s in summary),
        worsened_arms=sum(s['worsened'] for s in summary),
        unavailable_arms=sum(s['unavailable'] for s in summary))
    save(P/'summary.json', dict(totals=totals, views=summary,
        independent_confirmation=False, physical_source_truth=False))
    text = ['# Focused reconciliation sensitivity — 23 October 2026', '',
        'All eight acquisitions / nine views and all 81 specified arms complete. '
        'This is a post-hoc sensitivity of exposed observations, not independent confirmation. '
        'The original packages and failed outcomes are unchanged.', '',
        '## Does restricting the code statistic to the reference population narrow the gap?', '',
        'The comparison below uses centered float64 for both the full and selected code populations. '
        'It separates source selection from historical phase-precision differences. '
        'Within each setting, selected code and diagonal use the same localization rows and Fourier '
        'terms. Full and selected populations can target different weighted attenuations; '
        'the fitted equivalent-T gap is descriptive. '
        'Empty, nonpositive pooled-ratio and optimizer-boundary comparisons are unavailable in '
        'the counts below, with raw fit diagnostics retained. No equivalence or physical '
        'accuracy follows from a smaller gap.', '',
        '| View | Narrow / worsen / unavailable | Sources across grid | Selected code-input rows | Original-row Jaccard | Matched fit gap range (nm²) |',
        '|---|---:|---:|---:|---:|---:|']
    for s in summary:
        rr = s['ranges']; coverage = rr['localization_coverage']; overlap = rr['original_row_jaccard']
        gap = rr['selected_to_diagonal_fit_gap']
        text.append(f"| {s['tag']} | {s['narrowed']}/{s['worsened']}/{s['unavailable']} | "
            f"{rr['sources'][0]}–{rr['sources'][1]} | {100*coverage[0]:.1f}–{100*coverage[1]:.1f}% | "
            f"{overlap[0]:.3f}–{overlap[1]:.3f} | {gap[0]:.3f}–{gap[1]:.3f} |")
    text.extend(['', f"Gap narrowing occurs in {totals['narrowed_arms']}/81 arms and all nine settings "
        f"of {totals['all_settings_narrow']}/9 views. All nine settings worsen the gap in "
        f"{totals['all_settings_worsen']}/9 views; {totals['unavailable_arms']} arms are unavailable. "
        'These correlated grid cells are not independent replicates or hypothesis tests.', '',
        'Optimizer-boundary counts by view: '+', '.join(f"{s['tag']} {s['arithmetic_boundary_arms']}/9" for s in summary)+'. '
        'Both frequencies have usable selected-code and diagonal Fieller sets in '+', '.join(
            f"{s['tag']} {s['both_frequencies_fieller_usable_arms']}/9" for s in summary)+'. '
        'Point-fit gap narrowing alone is not supported precision.', '',
        '## Residual differences at the actual frequencies', '',
        'Counts below are conditional paired pointwise 95% bootstrap intervals excluding zero, '
        'not multiplicity-adjusted findings. Zero-spanning sets do not establish equivalence. '
        'Intervals are withheld if fewer than 95% of 2,000 draws have positive denominators/ratios.', '',
        '| View | q=.10 positive / negative / unavailable | q=.15 positive / negative / unavailable |',
        '|---|---:|---:|'])
    for s in summary:
        cells = [' / '.join(str(v[k]) for k in ('positive_pointwise_intervals',
            'negative_pointwise_intervals', 'unavailable_intervals')) for v in s['residuals']]
        text.append(f"| {s['tag']} | {cells[0]} | {cells[1]} |")
    text.extend(['', '## Quadratic weighting sensitivity', '',
        'Ranges cover all nine settings. These are different weighted finite-panel '
        'increment products, not interchangeable physical variances. No weights are '
        'chosen to favor agreement with a Fourier fit.', '',
        '| View | Uniform range | Event-count range | Localization-count range |',
        '|---|---:|---:|---:|'])
    for view in views:
        cells = []
        for name in WEIGHTS:
            values = [a['quadratic'][name]['point'] for a in view['arms']
                if a['quadratic'][name]['point'] is not None]
            cells.append(f'[{min(values):.3f}, {max(values):.3f}]' if values else 'unavailable')
        text.append(f"| {view['tag']} | {cells[0]} | {cells[1]} | {cells[2]} |")
    singleton = json.loads((P/'singleton-diagnostic.json').read_text())
    text.extend(['', '## Algebraic agreement diagnostic — post hoc', '',
        'After the primary grid completed, very sparse SDC selections motivated an '
        'additional diagnostic. With exactly one inferred source in a coarse group, '
        'the code and assigned diagonal are identical algebraically; agreement there '
        'cannot validate physical assumptions. An analytic control covers every code '
        'and unequal panel counts. The identity also holds for every retained singleton '
        'group within numerical tolerance.', '',
        'Restrict all three contrasts (full code, selected code and diagonal) to the '
        'same coarse groups containing more than one selected inferred source. This '
        'changes the conditional population and is secondary, not a replacement for '
        'the specified 81-arm result. Retain empty or unsupported cases.', '',
        '| View | Multi-source-only narrow / worsen / unavailable |', '|---|---:|'])
    for view in singleton['views']:
        counts = [sum(a['gap_narrows'] is value for a in view['arms']) for value in (True, False, None)]
        text.append(f"| {view['tag']} | {counts[0]}/{counts[1]}/{counts[2]} |")
    st = singleton['totals']
    text.extend(['', f"This diagnostic narrows {st['narrowed']}/81, worsens {st['worsened']}/81 "
        f"and leaves {st['unavailable']}/81 unsupported: {st['narrowed']}/{81-st['unavailable']} "
        'supported comparisons narrow. The primary tendency therefore is not explained '
        'solely by deterministic singleton agreement. It still does not establish '
        'physical assignments, equivalence or an all-source accuracy certificate.'])
    text.extend(['', '## Conditional reference-response uncertainty', '',
        'The original fitted responses replay for all six CD20 preparations. '
        'Observation and response-coefficient intervals are propagated by the specified conservative '
        'Bonferroni rectangle. Each noise law remains separate; their envelope is a sensitivity '
        'envelope, not a confidence interval for every physical model. '
        'Calibrated responses apply only to the original QC. NPC and SDC lack these matched '
        'response controls and are explicitly unsupported for this transformation.', '',
        '| Preparation | Weight | Original qualified | Gaussian propagated interval | Residual propagated interval |',
        '|---|---|---:|---|---|'])
    for r in responses['results']:
        for name, w in r.get('weights', {}).items():
            cells = []
            for model in ('gaussian', 'residual'):
                ci = w['models'][model]['propagated']['interval']
                cells.append(f'[{ci[0]:.3f}, {ci[1]:.3f}]' if ci is not None else 'unsupported/unbounded')
            text.append(f"| {r['tag']} | {name} | {w['original_response_qualified']} | {cells[0]} | {cells[1]} |")
    text.extend(['', '## Interpretation limits', '',
        'The primary source of uncertainty tested here is dependence on inferred clustering and '
        'eligibility, plus finite observation and response-fit uncertainty. Binding-event formation, '
        'precision filtering and temporal panels are fixed. There is no independent assignment truth '
        'and no validation of unselected sources. Weight changes alter the quadratic estimand; they '
        'are not competing estimates of one universal variance. Spatial blocks do not supply '
        'biological replication. The response intervals remain conditional on inferred templates, '
        'original QC and the two specified noise laws.', '',
        'All per-group contributions, complete Fieller set topologies, leave-one-block-out results, '
        'invalid bootstrap counts, baseline replays and raw-input hashes are retained. '
        'See PROTOCOL.md, all-arms.csv, summary.json, response-uncertainty.json and verification.json.'])
    (P/'RESULTS.md').write_text('\n'.join(text)+'\n')
    print(json.dumps(totals, indent=2))


if __name__ == '__main__':
    run()
