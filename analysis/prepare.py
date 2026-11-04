"""Prepare lossless figure tables and a source map from retained evidence only."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TAGS = ('npc', 'edf10rep1', 'edf10rep2', 'fig4untreated', 'rtx_rep1',
        'rtx_rep2', 'rtx_rep3', 'sdc_corrected', 'sdc_raw')


def read(name):
    return json.loads((ROOT / name).read_text())


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def prepare(out):
    out.mkdir(parents=True, exist_ok=False)
    sources = ['outputs/reference-reconciliation-oct07/comparison.json',
        'outputs/reconciliation-sensitivity-oct23/summary.json',
        'outputs/reconciliation-sensitivity-oct23/singleton-diagnostic.json',
        'methods/temporal-variants.md',
        'research/finite-code-covariance/THEORY.md',
        'research/finite-code-covariance/report-final/summary.json',
        'research/finite-code-covariance/report-final/all-comparisons.csv']
    original = read('outputs/reference-reconciliation-oct07/comparison.json')
    assert tuple(r['tag'] for r in original['acquisitions']) == TAGS
    views = []
    for r in original['acquisitions']:
        full, selected, diagonal = [r[k]['total_arithmetic'] for k in
                                   ('full_code', 'selected_code', 'assigned_diagonal')]
        views.append(dict(view=r['tag'], full=full, selected=selected, diagonal=diagonal,
                          signed_full_gap=full-diagonal, signed_selected_gap=selected-diagonal,
                          absolute_full_gap=abs(full-diagonal),
                          absolute_selected_gap=abs(selected-diagonal),
                          code_input_coverage=r['localization_coverage'],
                          frequency_differences=r['spatial_intervals']['paired_selected_difference']))
    write(out/'figure2-all-views.json', dict(
        units='nm2-equivalent descriptive point fits; not physical accuracy',
        interval_scope='central95% of sign-valid paired spatial-block bootstrap draws; conditional on inferred selection and validity; physical95% coverage unestablished',
        distinct_acquisitions=8, views=views))
    rows = []
    for tag in TAGS:
        name = f'outputs/reconciliation-sensitivity-oct23/result-{tag}.json'
        sources.append(name)
        for a in read(name)['arms']:
            rows.append(dict(view=tag, setting=a['name'], sources=a['sources'],
                code_input_coverage=a['localization_coverage'],
                full_gap=a['full_to_diagonal_fit_gap'],
                matched_gap=a['selected_to_diagonal_fit_gap'], narrows=a['fit_gap_narrows'],
                both_frequencies_usable=all(f['usable'] for key in ('selected_code', 'diagonal')
                    for f in a[key]['uncertainty']['per_frequency']),
                original_row_jaccard=a['original_row_jaccard']))
    with (out/'figure3-all81-settings.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    write(out/'figure3-summary.json', dict(
        original=read('outputs/reconciliation-sensitivity-oct23/summary.json'),
        singleton_diagnostic=read('outputs/reconciliation-sensitivity-oct23/singleton-diagnostic.json')['totals']))
    covariance = read('research/finite-code-covariance/report-final/summary.json')
    ratios = []
    for path in sorted((ROOT/'research/finite-code-covariance/results-attempt-01').glob('*.json')):
        d = json.loads(path.read_text())
        sources.append(str(path.relative_to(ROOT)))
        if d.get('theorem_applicable'):
            v = d['verification']
            for i, (m, exact) in enumerate(zip(v['mc_second_moment'], v['predicted_covariance'])):
                ratios.append(dict(world=path.name, component=i, mc_variance=m[i],
                                   exact_variance=exact[i], ratio=m[i]/exact[i]))
    write(out/'figure4-complete.json', dict(summary=covariance, component_variance_ratios=ratios))
    # Transfer the complete long condition table; nothing is selected or re-fit.
    temporal=(ROOT/'methods/temporal-variants.md').read_text().split(
        '## Complete observed-count coded primary table',1)[1]
    table=[]
    for line in temporal.splitlines():
        if not line.startswith('|'):continue
        cells=[c.strip() for c in line.strip().strip('|').split('|')]
        if all(set(c)<=set(':- ') for c in cells):continue
        table.append(cells)
    assert table and all(len(r)==len(table[0]) for r in table)
    with (out/'supplementary-data1-temporal-conditions.csv').open('w',newline='') as stream:
        csv.writer(stream,lineterminator='\n').writerows(table)
    assert len(rows) == 81 and len(ratios) == 864
    facts = dict(views=len(views), distinct_acquisitions=8,
        original_narrowed=sum(v['absolute_selected_gap'] < v['absolute_full_gap'] for v in views),
        sensitivity_narrowed=sum(r['narrows'] for r in rows),
        sensitivity_widened=sum(not r['narrows'] for r in rows),
        both_frequencies_usable=sum(r['both_frequencies_usable'] for r in rows),
        selected_coverage_min=min(r['code_input_coverage'] for r in rows),
        selected_coverage_max=max(r['code_input_coverage'] for r in rows),
        singleton=read('outputs/reconciliation-sensitivity-oct23/singleton-diagnostic.json')['totals'],
        covariance_worlds=covariance['worlds'], covariance_applicable=covariance['exogenous_worlds'],
        covariance_excluded=covariance['excluded_code_offset_worlds'],
        finite_frequency_variant_cases=covariance['finite_frequency_variant_cases'])
    write(out/'facts.json', facts)
    manifest = [dict(path=p, bytes=(ROOT/p).stat().st_size,
                     sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()) for p in sorted(set(sources))]
    write(out/'source-pins.json', dict(files=manifest, new_fitting=False, new_outcome_study=False))
    print(json.dumps(facts, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path, help='New directory only')
    prepare(parser.parse_args().output.resolve())
