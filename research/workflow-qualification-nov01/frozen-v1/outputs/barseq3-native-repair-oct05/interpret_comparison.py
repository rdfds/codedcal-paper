"""A claim-level synthesis that retains successful alternatives and transfer failures."""
import numpy as np
from condition import HERE,ROOT,read,sha


def main():
    primary=read(HERE/'recovery/summary.json');inner=read(HERE/'fit-support/summary.json')
    assert primary['full_cohort_comparison_supported']
    assert inner['counts']['missing']==0
    native_check=read(HERE/'verification.json');inner_check=read(HERE/'fit-support/verification.json')
    global_check=read(HERE/'assignment-verification.json')
    assert native_check['recovery_arms_verified']==60
    assert inner_check['completed_verified']+inner_check['failures_retained']==60
    lines=['# Earned public-data conclusion — 6 October 2026','',
        '**The added evidence supports a bounded comparative conclusion, not a validated automatic recommendation or a confident Nature Methods acceptance prediction.**','',
        'At identical complete injected offset distributions, assigning error by round versus color changes which tested calibration workflow restores the original molecular outputs. Round-only correction restores both common cases exactly; it does not restore the color cases. The tested barcode-aware author optimizer is far stronger than image-template alignment on one field, but that success does not transfer to the other. Low point blank proxies do not certify consistency. These statements retain every strong-alternative result and do not require our estimator to win.','',
        'The external comparison is on two fields of one BARseq3 preparation. The preceding annotated CSHL benchmark is a different public study, with a different and incomplete reference endpoint. Their combination supports conditional calibration/decoder dependence across public releases; it does not establish a universal quantitative recovery law across independent animals, prospective physical validation, or a bridge from the temporal contrast T to RNA accuracy.','',
        '## Subsequent exposed-data repair','',
        'The original failures below remain valid outcomes of their frozen workflows. Subsequent objective diagnosis finds a better oracle-initialized basin; pooling alone still fails field one. A structured global group search, using observed images and the codebook without injected offsets, then restores byte-identical native inputs and calls in all twelve field/case/decoder combinations (F1 1.000, composition TV 0). All twenty-two synthetic controls pass. Optional author refinement harms one common-round case and is excluded from the selected future policy. Source Data 66 and Supplementary Figure 6 retain the complete diagnostic and repair ledger; REPAIR_RESULTS.md gives its bounded scope. This repairs an exposed computational failure, not the original selector or independent transfer.','',
        'Two previously unexposed D076 preparations are metadata-qualified candidates for transfer testing. D078 is excluded because the project already used its raw release in Source Data 23. The public image server times out; no D076 pixels have been acquired or examined. Only the calibration core is frozen; the source adapter and complete independent workflow are not yet qualified.','',
        '## Complete matched-distribution comparison','',
        '| Field | Decoder | Common-case strong round F1 | Color-case strong round F1 | Color-case template F1 | Color-case author model F1 | Color-case oracle F1 |',
        '|---|---|---:|---:|---:|---:|---:|']
    for i in [0,1]:
        for d in ['correlation','regression']:
            cases=primary['fields'][str(i)]['decoders'][d]['cases'];common=cases[1]['scores'];color=cases[2]['scores']
            lines.append(f"| {i} | {d} | {common['common_iterative']['5']:.4f} | {color['common_iterative']['5']:.4f} | {color['population_template']['5']:.4f} | {color['barcode_model']['5']:.4f} | {color['oracle']['5']:.4f} |")
    lines+=['','All 60/60 arms complete without execution failures; quality misses remain failures of their specific targets. Decoder references use the same window as every corrected arm. Native correlation calls independently reproduce in NumPy. All integer inverse controls reproduce exact call lists, not just similar summary statistics.','',
        '## Border and context sensitivity','',
        f"All sixty planned follow-up records are retained: {inner['counts']['completed']} completed, {inner['counts']['failed']} failed, {inner['counts']['missing']} missing. This was specified after primary development exposure, before follow-up outcomes; it is post hoc sensitivity rather than another confirmation cohort.",'',
        '| Field | Decoder | Original color template F1 | Interior color template F1 | Original color model F1 | Interior color model F1 | Original chosen F1 | Interior chosen F1 |',
        '|---|---|---:|---:|---:|---:|---:|---:|']
    for i in [0,1]:
        for d in ['correlation','regression']:
            before=primary['fields'][str(i)]['decoders'][d]['cases'][2]
            after=inner['fields'][str(i)][d]['cases'][2]
            if len(after['scores'])<5:lines.append(f'| {i} | {d} | — | Failed fit retained | — | — | — | — |');continue
            a=before['scores'];b=after['scores']
            lines.append(f"| {i} | {d} | {a['population_template']['5']:.4f} | {b['population_template']['5']:.4f} | {a['barcode_model']['5']:.4f} | {b['barcode_model']['5']:.4f} | {a[before['recommend']]['5']:.4f} | {b[after['recommend']]['5']:.4f} |")
    lines+=['','Interior fits exclude artificial empty borders and check observed support at every template iteration. The smaller domain also changes useful information, so the difference cannot be attributed exclusively to padding. All zero/common cases, unnecessary corrections, three radii and complete method ledgers remain. The original cohort is not overwritten.','',
        '## What survives scrutiny','',
        '- Readout allocation matters at exactly the same full injected displacement distribution, in both external fields and both actual decoder families. This concerns imposed errors, not equal post-correction residuals or equal measured temporal T.',
        '- Output-count and blank gates can pass an arm with almost no native identity consistency. On external field one, strong round-only regression returns 437 nonblank calls, zero blanks, and F1 0.00795 against 66 native reference calls. An equally optimal spatial assignment gives F1 0.01193; gene-label-constrained matching at five pixels gives 0.13519, while composition TV is 0.87643. The severe output change persists, but spatial ambiguity limits the precise identity numerator. This is a controlled invariance failure, not a measured real-data RNA false-discovery rate.',
        '- The restricted-support correlation bound and selected actual-kernel absent-class witness make specific codebook/decoder mechanisms testable. The real data do not justify a universal isolated-source model: overlap, PSF tails and background remain.',
        '- The existing early MERFISH calibration/composition evidence remains useful, with strong native/image/landmark alternatives retained. This comparative extension does not supply unique early-method superiority.','',
        '## What has not been earned','',
        'A safe workflow selector, a universally accurate calibration certificate, an independently confirmed physical temporal intervention, and a consistent Nature Methods acceptance case are not established. Registration and barcode-dependent errors were already known; the exact matched-distribution comparison adds a narrower quantified boundary. Its importance beyond those established observations remains an editorial and scientific risk. More reused perturbations or favorable threshold tuning cannot turn exposed data into independent confirmation.','',
        'The four reserved kidney fields remain unopened. Failed S3 and ISS tests, two temporal-estimator modality classes, dense/sparse fiducial qualifications and the simulation-only randomized MERFISH retrofit remain unchanged. No private data, outside review, outreach, paid run, merge or submission was used for this extension.','',
        '## Verification and reproducibility','',
        f"Retained-call checks verify 120 integer-control arms and 60 external primary arms at three radii. The interior check verifies {inner_check['completed_verified']} completed records, {inner_check['failures_retained']} failures and {inner_check['exact_input_reuses_verified']} byte-identical input reuses. Whole-domain matching verifies cardinality and minimum distance for {global_check['completed_arms']} arms, retaining {global_check['equal_cost_identity_differences']} equal-cost identity differences. The largest observed F1 change between equally optimal spatial matchings is {global_check['largest_observed_equal_cost_f1_change']:.5f}; this compares two assignments and is not a bound over every tied solution.",'',
        'Source Data 65 indexes hashes, frozen protocols, environments, runtime components, support/blank/output misses, counterexample limits and public archive audits. All six main figures remain; Supplementary Figure 5 displays every primary condition. This report SHA is generated from its committed analysis script and retained source summaries.', '',
        f"Analysis script SHA256: `{sha(__file__)}`."]
    (HERE/'CONCLUSION.md').write_text('\n'.join(lines)+'\n')
    print('Wrote complete claim-level comparison synthesis',flush=True)

if __name__=='__main__':main()
