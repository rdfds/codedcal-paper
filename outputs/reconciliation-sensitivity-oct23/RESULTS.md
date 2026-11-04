# Focused reconciliation sensitivity — 23 October 2026

All eight acquisitions / nine views and all 81 specified arms complete. This is a post-hoc sensitivity of exposed observations, not independent confirmation. The original packages and failed outcomes are unchanged.

## Does restricting the code statistic to the reference population narrow the gap?

The comparison below uses centered float64 for both the full and selected code populations. It separates source selection from historical phase-precision differences. Within each setting, selected code and diagonal use the same localization rows and Fourier terms. Full and selected populations can target different weighted attenuations; the fitted equivalent-T gap is descriptive. Empty, nonpositive pooled-ratio and optimizer-boundary comparisons are unavailable in the counts below, with raw fit diagnostics retained. No equivalence or physical accuracy follows from a smaller gap.

| View | Narrow / worsen / unavailable | Sources across grid | Selected code-input rows | Original-row Jaccard | Matched fit gap range (nm²) |
|---|---:|---:|---:|---:|---:|
| npc | 4/5/0 | 387–4695 | 4.4–37.5% | 0.256–1.000 | 0.154–1.082 |
| edf10rep1 | 9/0/0 | 493–4205 | 4.4–28.8% | 0.318–1.000 | 0.017–0.622 |
| edf10rep2 | 9/0/0 | 2314–5983 | 28.6–61.1% | 0.620–1.000 | 0.726–1.920 |
| fig4untreated | 4/5/0 | 403–1663 | 12.6–49.4% | 0.396–1.000 | 0.017–2.014 |
| rtx_rep1 | 9/0/0 | 805–6045 | 7.2–45.6% | 0.304–1.000 | 0.018–0.472 |
| rtx_rep2 | 8/1/0 | 1434–4913 | 18.8–56.7% | 0.512–1.000 | 0.010–0.240 |
| rtx_rep3 | 9/0/0 | 1242–8134 | 8.6–47.5% | 0.329–1.000 | 0.028–0.238 |
| sdc_corrected | 7/2/0 | 15–2915 | 0.1–15.8% | 0.030–1.000 | 0.000–0.760 |
| sdc_raw | 9/0/0 | 11–2494 | 0.1–13.4% | 0.027–1.000 | 0.000–0.704 |

Gap narrowing occurs in 68/81 arms and all nine settings of 5/9 views. All nine settings worsen the gap in 0/9 views; 0 arms are unavailable. These correlated grid cells are not independent replicates or hypothesis tests.

Optimizer-boundary counts by view: npc 0/9, edf10rep1 0/9, edf10rep2 0/9, fig4untreated 0/9, rtx_rep1 0/9, rtx_rep2 0/9, rtx_rep3 0/9, sdc_corrected 0/9, sdc_raw 0/9. Both frequencies have usable selected-code and diagonal Fieller sets in npc 9/9, edf10rep1 9/9, edf10rep2 9/9, fig4untreated 9/9, rtx_rep1 9/9, rtx_rep2 9/9, rtx_rep3 9/9, sdc_corrected 8/9, sdc_raw 7/9. Point-fit gap narrowing alone is not supported precision.

## Residual differences at the actual frequencies

Counts below are conditional paired pointwise 95% bootstrap intervals excluding zero, not multiplicity-adjusted findings. Zero-spanning sets do not establish equivalence. Intervals are withheld if fewer than 95% of 2,000 draws have positive denominators/ratios.

| View | q=.10 positive / negative / unavailable | q=.15 positive / negative / unavailable |
|---|---:|---:|
| npc | 0 / 0 / 0 | 5 / 0 / 0 |
| edf10rep1 | 0 / 0 / 0 | 0 / 0 / 0 |
| edf10rep2 | 4 / 0 / 0 | 9 / 0 / 0 |
| fig4untreated | 1 / 0 / 0 | 0 / 0 / 0 |
| rtx_rep1 | 1 / 0 / 0 | 0 / 0 / 0 |
| rtx_rep2 | 0 / 0 / 0 | 0 / 0 / 0 |
| rtx_rep3 | 0 / 0 / 0 | 0 / 0 / 0 |
| sdc_corrected | 1 / 0 / 0 | 1 / 0 / 0 |
| sdc_raw | 0 / 1 / 0 | 0 / 0 / 0 |

## Quadratic weighting sensitivity

Ranges cover all nine settings. These are different weighted finite-panel increment products, not interchangeable physical variances. No weights are chosen to favor agreement with a Fourier fit.

| View | Uniform range | Event-count range | Localization-count range |
|---|---:|---:|---:|
| npc | [0.593, 1.281] | [0.713, 1.217] | [0.669, 1.197] |
| edf10rep1 | [1.529, 1.889] | [1.642, 2.291] | [1.539, 2.455] |
| edf10rep2 | [1.663, 1.913] | [1.344, 1.973] | [1.986, 2.610] |
| fig4untreated | [3.127, 3.662] | [2.810, 3.818] | [2.819, 3.887] |
| rtx_rep1 | [2.065, 3.502] | [2.631, 4.555] | [3.095, 4.958] |
| rtx_rep2 | [1.101, 1.651] | [0.926, 1.928] | [0.949, 2.186] |
| rtx_rep3 | [1.123, 1.672] | [1.117, 1.621] | [1.232, 1.850] |
| sdc_corrected | [2.339, 6.048] | [1.072, 6.125] | [0.970, 9.002] |
| sdc_raw | [4.655, 17.284] | [4.272, 15.658] | [4.670, 14.918] |

## Algebraic agreement diagnostic — post hoc

After the primary grid completed, very sparse SDC selections motivated an additional diagnostic. With exactly one inferred source in a coarse group, the code and assigned diagonal are identical algebraically; agreement there cannot validate physical assumptions. An analytic control covers every code and unequal panel counts. The identity also holds for every retained singleton group within numerical tolerance.

Restrict all three contrasts (full code, selected code and diagonal) to the same coarse groups containing more than one selected inferred source. This changes the conditional population and is secondary, not a replacement for the specified 81-arm result. Retain empty or unsupported cases.

| View | Multi-source-only narrow / worsen / unavailable |
|---|---:|
| npc | 4/5/0 |
| edf10rep1 | 9/0/0 |
| edf10rep2 | 9/0/0 |
| fig4untreated | 5/4/0 |
| rtx_rep1 | 9/0/0 |
| rtx_rep2 | 8/1/0 |
| rtx_rep3 | 9/0/0 |
| sdc_corrected | 5/2/2 |
| sdc_raw | 7/0/2 |

This diagnostic narrows 65/81, worsens 12/81 and leaves 4/81 unsupported: 65/77 supported comparisons narrow. The primary tendency therefore is not explained solely by deterministic singleton agreement. It still does not establish physical assignments, equivalence or an all-source accuracy certificate.

## Conditional reference-response uncertainty

The original fitted responses replay for all six CD20 preparations. Observation and response-coefficient intervals are propagated by the specified conservative Bonferroni rectangle. Each noise law remains separate; their envelope is a sensitivity envelope, not a confidence interval for every physical model. Calibrated responses apply only to the original QC. NPC and SDC lack these matched response controls and are explicitly unsupported for this transformation.

| Preparation | Weight | Original qualified | Gaussian propagated interval | Residual propagated interval |
|---|---|---:|---|---|
| edf10rep1 | unweighted | True | [0.735, 9.644] | [0.468, 9.516] |
| edf10rep1 | event_count | True | [1.231, 8.565] | [1.212, 8.978] |
| edf10rep1 | localization_count | True | [1.056, 35.084] | [1.094, 37.275] |
| edf10rep2 | unweighted | True | [0.622, 2.307] | [0.645, 2.382] |
| edf10rep2 | event_count | True | [1.064, 2.547] | [0.916, 2.366] |
| edf10rep2 | localization_count | True | [1.049, 3.734] | [0.876, 3.711] |
| fig4untreated | unweighted | True | [2.226, 4.921] | [2.478, 5.251] |
| fig4untreated | event_count | True | [1.781, 4.778] | [1.905, 5.151] |
| fig4untreated | localization_count | False | [1.792, 5.019] | [2.058, 5.629] |
| rtx_rep1 | unweighted | True | [1.687, 18.035] | [1.567, 18.449] |
| rtx_rep1 | event_count | False | [2.094, 8.925] | [2.071, 9.595] |
| rtx_rep1 | localization_count | False | [2.661, 13.440] | [2.834, 15.569] |
| rtx_rep2 | unweighted | True | [0.896, 1.795] | [0.797, 1.770] |
| rtx_rep2 | event_count | True | [0.901, 3.961] | [0.861, 4.392] |
| rtx_rep2 | localization_count | False | [1.029, 4.246] | [0.968, 4.276] |
| rtx_rep3 | unweighted | True | [0.326, 1.996] | [0.366, 2.247] |
| rtx_rep3 | event_count | False | [0.477, 2.160] | [0.425, 2.446] |
| rtx_rep3 | localization_count | False | [0.531, 2.867] | [0.544, 3.231] |

## Interpretation limits

The primary source of uncertainty tested here is dependence on inferred clustering and eligibility, plus finite observation and response-fit uncertainty. Binding-event formation, precision filtering and temporal panels are fixed. There is no independent assignment truth and no validation of unselected sources. Weight changes alter the quadratic estimand; they are not competing estimates of one universal variance. Spatial blocks do not supply biological replication. The response intervals remain conditional on inferred templates, original QC and the two specified noise laws.

All per-group contributions, complete Fieller set topologies, leave-one-block-out results, invalid bootstrap counts, baseline replays and raw-input hashes are retained. See PROTOCOL.md, all-arms.csv, summary.json, response-uncertainty.json and verification.json.
