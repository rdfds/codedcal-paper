"""Summarize every frozen world; no selection of favorable comparisons."""
from pathlib import Path
import argparse,json,csv
from collections import defaultdict
import numpy as np


def main():
    p=argparse.ArgumentParser();p.add_argument('--prefix',required=True);p.add_argument('--finite',required=True);p.add_argument('--output',required=True);args=p.parse_args()
    prefix,finite,out=Path(args.prefix),Path(args.finite),Path(args.output)
    if out.exists():raise ValueError('new output directory required')
    assert json.loads((prefix/'state.json').read_text())['stage']=='complete'
    assert json.loads((finite/'state.json').read_text())['stage']=='complete'
    out.mkdir(parents=True)
    covs=[json.loads(f.read_text()) for f in sorted(prefix.glob('*.json')) if f.name!='state.json']
    budgets=[json.loads(f.read_text()) for f in sorted(finite.glob('*.json')) if f.name!='state.json']
    assert len(covs)==126 and len(budgets)==108
    rows=[];all_rows=[]
    for d in budgets:
        for r in d['comparisons']:
            row={k:v for k,v in r.items() if k not in ('cases',)}
            row.update(name=d['name'],tag=d['name'].rsplit('-',2)[0],model=d['name'].rsplit('-',2)[1],world=d['name'].rsplit('-',2)[2])
            all_rows.append(row)
            if r['supported']:rows.append(row)
    keys=sorted(set().union(*(r.keys() for r in all_rows)))
    with (out/'all-comparisons.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(all_rows)
    index={(r['name'],r['q'],r['variant'],r['K']):r for r in rows}
    summaries=[]
    for variant in ('arithmetic','geometric'):
        for k in (4,8,16):
            selected=[r for r in rows if r['K']==k and r['variant']==variant]
            ratios=[];rmse=[];better=worse=tie=0
            for r in selected:
                ref=index[r['name'],r['q'],r['variant'],4]
                ratios.append(r['predicted_influence_variance']/ref['predicted_influence_variance'])
                if r['finite_available_rmse'] is not None and ref['finite_available_rmse'] not in (None,0):rmse.append(r['finite_available_rmse']/ref['finite_available_rmse'])
                better+=r['usable_covers']>ref['usable_covers'];worse+=r['usable_covers']<ref['usable_covers'];tie+=r['usable_covers']==ref['usable_covers']
            summaries.append(dict(variant=variant,K=k,supported_conditions=len(selected),unsupported_conditions=216-len(selected),
                median_predicted_variance_factor=float(np.median(ratios)),variance_improvements=sum(x<1-1e-12 for x in ratios),
                median_available_rmse_factor=float(np.median(rmse)),
                available=sum(r['available'] for r in selected),usable=sum(r['usable'] for r in selected),
                usable_covers=sum(r['usable_covers'] for r in selected),total=216*200,
                usable_coverage_better=better,usable_coverage_worse=worse,usable_coverage_tied=tie))
    verification=[d['verification'] for d in covs if d['theorem_applicable']]
    z=np.array([x for r in verification for x in r['variance_z']])
    data=dict(worlds=126,exogenous_worlds=108,excluded_code_offset_worlds=18,
              exact_prefix_known=2016,exact_prefix_coded=4032,
              finite_experiments=108*200*3,finite_frequency_variant_cases=108*200*3*4,
              supported_cases=sum(200 for r in all_rows if r['supported']),
              covariance_relative_frobenius_median=float(np.median([v['relative_frobenius_error'] for v in verification])),
              covariance_relative_frobenius_max=float(max(v['relative_frobenius_error'] for v in verification)),
              variance_z_abs_max=float(max(abs(z))),variance_z_abs_over_three=int(sum(abs(z)>3)),
              summary=summaries)
    (out/'summary.json').write_text(json.dumps(data,indent=2)+'\n')
    lines=['# Finite-code covariance and design extension — results', '',
           '27 October 2026. Bounded conditional theory/computation; original manuscript,',
           'Source Data 69–72, failed results and shipped estimator remain unchanged.', '',
           '## Main finding', '',
           'The exact covariance of the implemented source-coded Fourier moments separates',
           'into pair and triple Gram matrices with alphabet factors 1/(K−1) and',
           '1/((K−1)(K−2)). Physical-field variation is a separate floor. For a fixed',
           'scalar ratio influence, eight codes beat four at B=K×G exactly when',
           '**A < P/21 + 5Q/42**. These quantities require source contributions or a',
           'justified source model; this is not an observable source-blind recommendation.', '',
           'There is a real **model-conditional** variance/design improvement, but no',
           'universal improvement in usable measurement. Fewer independent groups can',
           'make Fieller sets unusable despite lower first-order variance. No physical',
           'hardware cost, new random-readout experiment or independent molecular truth',
           'has been established.', '',
           '## Verification and complete scope', '',
           '- All 126 historical worlds: 2,016 known-source and 4,032 coded prefix',
           '  contribution vectors reproduced bit-for-bit.',
           '- All nine views/eight acquisitions and both noise laws included.',
           '- 108 exogenous worlds receive conditional covariance/design evaluation;',
           '  all 18 code-dependent-offset worlds remain explicit theorem exclusions.',
           '- 512 independent code assignments per fixed sixteen-group observation',
           '  prefix check covariance; these are numerical replicates, not acquisitions.',
           f'- Median relative matrix discrepancy is {data["covariance_relative_frobenius_median"]:.1%}; maximum {data["covariance_relative_frobenius_max"]:.1%}.',
           f'- {data["variance_z_abs_over_three"]}/864 component variance differences exceed three *estimated* Monte Carlo standard errors; maximum {data["variance_z_abs_max"]:.2f}.',
           '  These discrepancies are retained; estimated Monte Carlo errors are not',
           '  treated as calibrated Gaussian agreement tests. No extra draws were added',
           '  after disagreement. Exact proof/enumeration is a separate check.',
           '- 64,800 finite-budget simulations (200 per K/world), K=4,8,16 and',
           '  B=256, give 259,200 planned frequency/variant cases. No K is selected',
           '  from the outcomes. Every sign/support failure and confidence-set type',
           '  is retained. Original failures are not replaced by these simulations.', '',
           '## Equal-budget results', '',
           'The observation law is uniform on the SAME fixed sixteen reconstructed',
           'groups, with fresh independent permanent codes. Both terms of total',
           'covariance and each variant’s target use that same finite-state population.',
           'Counts, photons and kernels are held unchanged by stipulation as K changes.',
           'This establishes a conditional computation, not a measured acquisition design.', '',
           '| Variant | K / groups | Median predicted variance / K4 | Conditions with lower variance | Median available-ratio RMSE / K4 | Usable + covering / all planned calls | Conditions better / worse / tied usable coverage |',
           '|---|---:|---:|---:|---:|---:|---:|']
    for r in summaries:
        lines.append(f'| {r["variant"]} | {r["K"]} / {256//r["K"]} | {r["median_predicted_variance_factor"]:.3f} | {r["variance_improvements"]}/{r["supported_conditions"]} | {r["median_available_rmse_factor"]:.3f} | {r["usable_covers"]}/{r["total"]} ({r["usable_covers"]/r["total"]:.1%}) | {r["usable_coverage_better"]} / {r["usable_coverage_worse"]} / {r["usable_coverage_tied"]} |')
    lines+=['', 'Available-ratio RMSE is conditional on finite positive point support; it is',
            'not unconditional risk and must be read with all availability/usability',
            'counts in all-comparisons.csv. World/frequency conditions overlap and are',
            'not independent biological replications. No acceptance or coverage guarantee',
            'is inferred from pooled descriptive percentages.', '',
            '### All views: homogeneous world, arithmetic, both noise laws and frequencies', '',
            'Each entry aggregates four conditions / 800 simulated calls; model-specific',
            'details and geometric results remain in all-comparisons.csv. This table is',
            'a retrospective descriptive display; complete worlds remain in the CSV.', '',
            '| View | K8 median variance factor | K4 usable coverage | K8 usable coverage | K16 usable coverage |',
            '|---|---:|---:|---:|---:|']
    for tag in ('npc','edf10rep1','edf10rep2','fig4untreated','rtx_rep1','rtx_rep2','rtx_rep3','sdc_corrected','sdc_raw'):
        selected=[r for r in rows if r['tag']==tag and r['world']=='homogeneous' and r['variant']=='arithmetic']
        ratio=[r['predicted_influence_variance']/index[r['name'],r['q'],'arithmetic',4]['predicted_influence_variance'] for r in selected if r['K']==8]
        counts=[sum(r['usable_covers'] for r in selected if r['K']==k) for k in (4,8,16)]
        lines.append(f'| {tag} | {np.median(ratio):.3f} | {counts[0]}/800 | {counts[1]}/800 | {counts[2]}/800 |')
    lines+=['', '## What the result does and does not repair', '',
            '1. The joint-covariance calculation is now about the actual implemented',
            '   finite-phase moments, including unequal counts and dependence across',
            '   frequencies/windows/angles. It is not a denominator-noise argument.',
            '2. The restricted optimality proof now uses the correct five code-label',
            '   equality orbits. Arbitrary-contribution unbiasedness forces the three',
            '   repeated-position coefficients to agree. It is not a lower bound for',
            '   nonlinear estimators; the N=3 toy counterexample is retained.',
            '3. The generic mathematics is established prior art. The contraction is',
            '   exactly K times the classical third-cumulant polynomial across bins.',
            '   Dependent bins require the source-assignment proof, but do not make',
            '   the polynomial new. PRIOR_WORK.md limits the possible novelty claim.',
            '4. The first-pass population scores combine a sixteen-group code-variance',
            '   estimate with the full retained population’s physical floor. They are',
            '   preliminary descriptive estimates, not exact population optima. The',
            '   finite-budget verification deliberately uses a single consistent law.',
            '5. In particular the EDF10-rep2 prefix can miss rare dominant contributions.',
            '   Its new conditional performance cannot fix the original acquisition’s',
            '   poor Fieller calibration or erase the all-field failure. More codes',
            '   cannot remove a physical/contribution-weight floor.',
            '6. Neither covariance nor hypothetical code allocation validates physical',
            '   code exogeneity, labeling effects, source-specific kernel changes or',
            '   shared errors. Kernel-violation observable targets remain distinct from',
            '   physical attenuation. Randomized MERFISH remains simulation-only.', '',
            '## Decision', '',
            'Retain this as a bounded theory/design extension with a concrete conditional',
            'benefit and an explicit failure boundary. It improves explanation and',
            'measurement-design analysis; it does not establish a universal better',
            'estimator, a hardware-ready optimum or closure of the physical-validation',
            'gap. The inspected prior work does not justify marketing this as a new',
            'general cumulant or sketching principle. No manuscript claim upgrade or',
            'replacement of Fieller uncertainty is warranted automatically.', '',
            'Full derivation: ../THEORY.md. Prior-work comparison: ../PRIOR_WORK.md.',
            'Controls and independent recount: ../test_covariance.py and ../verify.py.',
            'Historical artifacts were read-only throughout.']
    (out/'RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(data))

if __name__=='__main__':main()
