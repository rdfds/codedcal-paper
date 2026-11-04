"""Render the complete comparison, including all model and acquisition failures."""
from common import P,ROOT,read

def main():
    report=read(str((P/'comparison.json').relative_to(ROOT)))
    a=report['acquisitions'];cd=[x for x in a if 'response_sensitivity' in x]
    lines=['# Existing coded-acquisition reference comparison — 7 October 2026','',
        'Eight distinct acquisitions / nine raw-or-corrected views are included. '
        'Figure4 RTX is the same acquisition as EDF11/EDF12 RTX replicate 1; '
        'all four picked-member hashes agree, and it counts once. '
        'All analyses are post hoc on exposed public data.','',
        '## What is now reconciled','',
        'The historical quadratic reference uses weighted binding-event means '
        'and event-count products on a clusterable subset. The code statistic '
        'uses unweighted localization Fourier sums and localization-count '
        'products on the original, broader grouped-code population. A matched assignment-assisted '
        'diagonal contrast and code replay on its exact retained localization '
        'rows separate these differences. Source selection remains '
        'oracle-dependent; this is a validation diagnostic, not a new '
        'assignment-free workflow or proof about unselected sources.','',
        'All eligible historical source rows, contrasts and count weights '
        'are replayed before interpreting new controls. The raw SDC view has '
        'no historical raw-oracle source file; its corrected oracle and both '
        'full-code moment arrays are audited. Legacy code fits below retain '
        'the historical bounded two-frequency fit solely for comparison; '
        'they are not unconditional physical variances.','',
        '## Interpretation','',
        'The reference comparison is now acquisition-complete and target-aware. '
        'For RTX replicate 1, the original grouped-code fit is 6.576 nm^2, '
        'while same-row code/diagonal fits are 3.550/3.531. Selection accounts '
        'for much of that old discrepancy. It does not improve every comparison: '
        'Figure4 untreated moves from an original fit of 2.541 to a same-row '
        'fit of 4.102, against a diagonal fit of 3.255.','',
        'EDF10 replicate 2 retains positive conditional code-minus-diagonal '
        'differences at both frequencies; corrected SDC retains one at q=0.10. '
        'Other same-row paired intervals span zero, which does not demonstrate '
        'equivalence. Oracle selection can affect code exogeneity, and inferred '
        'source assignments, heterogeneous motion and changing kernels remain '
        'possible limitations. These data cannot uniquely assign the remaining '
        'discrepancies to physical motion or estimator bias.','',
        'Declared conditional response qualification passes for 4/6 event-'
        'weighted and 2/6 localization-weighted CD20 references. An NPC-derived '
        'affine response is therefore not a generally validated transfer. '
        'The result repairs the validation comparison and its claims; it does '
        'not establish physical temporal truth, validation of unselected '
        'sources or a new independent benefit over strong alternatives.','',
        '## Every observed view','',
        '| View | Eligible sources | Code-input rows retained | Event-weight quadratic | Localization-weight quadratic | Original grouped-code fit | Same-row code fit | Same-row diagonal fit |',
        '|---|---:|---:|---:|---:|---:|---:|---:|']
    for x in a:
        v=x['estimates'];lines.append(f"| {x['tag']} | {x['passing_sources']:,} | {100*x['localization_coverage']:.1f}% | {v['event_count']['raw']:.3f} | {v['localization_count']['raw']:.3f} | {x['full_code']['total_arithmetic']:.3f} | {x['selected_code']['total_arithmetic']:.3f} | {x['assigned_diagonal']['total_arithmetic']:.3f} |")
    lines += ['', 'Coverage is relative to the original grouped-code input rows in valid '
        'temporal panels, not all localizations in the acquisition. '
        'Values have units nm^2; interpreting them as physical variances '
        'requires their respective observation models. '
        'Quadratic and finite-frequency fits remain different estimands. '
        'comparison.json retains raw estimates, both frequencies, Fieller '
        'ratio sets, paired spatial-block differences, invalid draws, source '
        'rejections, effective sample weights and every historical replay. '
        'A confidence interval spanning a difference of zero is not an '
        'equivalence result. Spatial blocks are not biological replicates.','',
        '## Same-row disagreements and uncertainty','',
        'The two-frequency bounded point fit can conceal unsupported ratios. '
        'Every frequency-specific difference and paired spatial-block interval '
        'is therefore reported separately. Intervals are conditional on the '
        'observed oracle selection. Intervals are pointwise and unadjusted for '
        'multiple comparisons; no familywise or equivalence claim was specified.','',
        '| View | q (nm^-1) | Code minus diagonal equivalent T | Conditional paired pointwise 95% interval | Valid bootstrap draws |',
        '|---|---:|---:|---|---:|']
    for x in a:
        for r in x['spatial_intervals']['paired_selected_difference']:
            point='unsupported' if r['code_minus_diagonal'] is None else f"{r['code_minus_diagonal']:.3f}"
            ci=r['conditional_interval'];interval='unsupported' if ci is None else f"[{ci[0]:.3f}, {ci[1]:.3f}]"
            lines.append(f"| {x['tag']} | {r['q']:.2f} | {point} | {interval} | {100*r['valid_fraction']:.1f}% |")
    lines += ['', 'Unsupported original grouped-code frequency-specific uncertainty:','']
    for x in a:
        for r in x['spatial_intervals']['full_code']['frequencies']:
            if r['T_interval'] is None:
                lines.append(f"- {x['tag']}, q={r['q']:.2f}: no supported finite equivalent-T interval; Fieller set `{r['kind']}`. The bounded two-frequency point fit is not an uncertainty certificate.")
    lines += ['', '## Count, noise and selection conditions','',
        'All candidate clusters precede the original source-QC gate. '
        'Quantiles below describe the retained empirical ingredients. '
        'Within-panel scatter can include merged emitters and motion; it is '
        'not independently known localization noise. Full per-panel counts, '
        'burst multiplicities, candidate-neighbour distances, frame boundaries '
        'and rejection reasons are retained in comparison.json and the NPZ tables.','',
        '| View | Candidate clusters | Median candidate → eligible events by panel | Median reported event variance (nm^2) | Median within-panel per-axis scatter (nm^2) | Event rows retained |',
        '|---|---:|---|---:|---:|---:|']
    for x in a:
        c=x['empirical_conditions'];pan='/'.join(f'{v:g}' for v in c['per_panel_event_count_quantiles'][1])
        pan+=' → '+'/'.join(f'{v:g}' for v in c['eligible_per_panel_event_count_quantiles'][1])
        lines.append(f"| {x['tag']} | {c['candidate_clusters']:,} | {pan} | {c['reported_event_variance_quantiles'][1]:.3f} | {c['within_panel_per_axis_scatter_quantiles'][1]:.3f} | {100*x['event_coverage']:.1f}% |")
    lines += ['', '## Preparation-specific reference response','',
        'Six CD20 acquisitions complete 96 control worlds and 384,000 group '
        'draws: T=0,1,3,5 nm^2; Gaussian and within-panel residual noise; '
        'two seeds; 4,000 group draws in every world. Inputs retain all '
        'candidate clusters before oracle source QC, their geometry, panel '
        'counts, noise and localization multiplicities. Frozen source QC '
        'is rerun after adding known window motion. Paired per-group '
        'influences preserve common random draws across T within each seed.','',
        '| Acquisition | Event-weight slopes (Gaussian / residual) | Event-weight response qualified | Localization-weight slopes (Gaussian / residual) | Localization-weight response qualified |',
        '|---|---:|---|---:|---|']
    for x in cd:
        event=x['response_sensitivity']['event_count'];loc=x['response_sensitivity']['localization_count']
        def slopes(v):return ' / '.join(f"{v['models'][m]['slope']:.3f}" for m in ('gaussian','residual'))
        lines.append(f"| {x['tag']} | {slopes(event)} | {event['qualified_conditional_response']} | {slopes(loc)} | {loc['qualified_conditional_response']} |")
    lines += ['', 'Failed criteria, including the unweighted reference sensitivity:','']
    for x in cd:
        for weight,r in x['response_sensitivity'].items():
            if r['qualified_conditional_response']:continue
            reasons=[]
            if r['relative_slope_difference']>.1:reasons.append(f"noise-model slope difference {100*r['relative_slope_difference']:.1f}% >10%")
            for model,f in r['models'].items():
                if f['slope_interval'][0]<=0:reasons.append(f'{model} slope interval includes zero')
                if f['max_affine_residual']>.2:reasons.append(f"{model} maximum affine residual {f['max_affine_residual']:.3f} >0.2 nm^2")
            lines.append(f"- {x['tag']}, {weight}: {'; '.join(reasons)}.")
    lines += ['', 'Qualification requires positive slope intervals, <=10% '
        'between-model slope difference, and <=0.2 nm^2 maximum affine '
        'residual. These are declared sensitivity criteria, not a field '
        'standard. Unqualified responses stay in the table and JSON. '
        'Intercepts, complete source-selection coverage, inferred-label '
        'impurity, model-specific transformed observations and uncertainty '
        'are retained. No model was selected to improve code/oracle agreement.','',
        'These are field-derived **conditional** controls. Inferred clusters '
        'may already merge real emitters; neither noise model establishes '
        'physical assignments, event independence or the actual displacement '
        'law. Transformed observation intervals do not include every '
        'reference-model uncertainty. A qualified affine response cannot '
        'certify the acquisition or turn quiet D into validation.','',
        '## Temporal support and the existing SDC calibration','',
        'NPC panels cover a 12,500-frame acquisition; CD20 panels are quarters '
        'of round-specific acquisitions ranging from 20,000 to 50,000 frames; '
        'SDC has three 4,200-frame panels and a final 4,250-frame panel, '
        'with 50-frame guard gaps over 17,000 '
        'frames. These are distinct finite temporal supports, not a universal '
        'common lag. Per-source mean event frames, exact localization and '
        'event panel counts, and round-specific boundaries are retained. '
        'No absolute synchronization of sequential rounds is inferred.','',
        'For noise-free panel means a,b,c,d, the quadratic reference equals '
        '[(a-c)(b-d)+(a-d)(b-c)]/4. It becomes a per-window variance only '
        'under the persistent independent-window model; general motion '
        'produces a signed increment-product functional. The finite-q '
        'contrast identifies attenuation under its separate assumptions.','',
        'The prior SDC simulation failed to reproduce every diagnostic and '
        'selection coverage. Its fitted intercept includes modeled settling '
        'motion, so subtracting it targets only the injected random component. '
        'It cannot be used as a fully matched reference for a code contrast '
        'retaining that physical motion. Both raw/corrected observations and '
        'the old conditional sensitivity remain. See ESTIMANDS.md.','',
        '## Numerical provenance and reproduction','',
        'The original CD20 code evaluated phases on global float32 positions, '
        'and the original SDC scalar phase operations also retained float32. '
        'New target-matching diagnostics use centered float64 positions. '
        'reconcile.py re-executes the historical code function to verify '
        'its original moments and records the separate numerical sensitivity '
        'without loosening replay tolerances. The first three control fields '
        'use dense legacy clustering; subsequent controls use a partition-'
        'equivalent sparse implementation, checked at the strict cutoff. '
        'The initial 400-draw pilot and an import failure before any controls '
        'remain; primary controls uniformly use 4,000 draws.','',
        'From the repository root, using ../merlin-venv.nosync/bin/python:', '',
        '```sh','python outputs/reference-reconciliation-oct07/fetch_existing.py',
        'python outputs/reference-reconciliation-oct07/pipeline.py',
        'python outputs/reference-reconciliation-oct07/test_controls.py',
        'python outputs/reference-reconciliation-oct07/reconcile.py',
        'python outputs/reference-reconciliation-oct07/write_report.py',
        'python outputs/reference-reconciliation-oct07/verify.py','```','',
        'NPC picked tables and SDC reduced projections are reused from the '
        'verified existing caches described in reference_characterize.py. '
        'The CD20 fetcher checks whole-member historical SHA-256, two '
        'extraction implementations and CRC, retaining one raw copy outside '
        'Git. Field-derived simulation ingredients, per-group controls and '
        'verification hashes are retained in this package.','']
    (P/'RESULTS.md').write_text('\n'.join(lines))

if __name__=='__main__':main()
