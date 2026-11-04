# Supplementary Table 2: existing coded-acquisition reference comparison

Source Data 69 includes every coded acquisition already used in the temporal
reference comparisons: eight distinct acquisitions, nine raw/corrected views.
Figure4 RTX and RTX replicate 1 have identical source-member hashes and count
once. These are post-hoc analyses of exposed data; no new preparation was
chosen for agreement. Historical source contrasts and original code moments
replay. All values and uncertainty retain their conditional model interpretation.
The code-input coverage denominator is the original grouped-code population,
not all localizations in the acquisition.

The exact estimand derivation is retained in
`../outputs/reference-reconciliation-oct07/ESTIMANDS.md`. Complete variants,
responses, rejected sources and input hashes are in comparison.json; independent
recounts are in verification.json. All interpretations and reproduction commands
appear in RESULTS.md.

## Every observed view

| View | Eligible sources | Code-input rows retained | Event-weight quadratic | Localization-weight quadratic | Original grouped-code fit | Same-row code fit | Same-row diagonal fit |
|---|---:|---:|---:|---:|---:|---:|---:|
| npc | 1,798 | 16.9% | 1.044 | 0.896 | 1.494 | 1.052 | 0.898 |
| edf10rep1 | 1,739 | 13.4% | 2.061 | 2.203 | 3.682 | 2.133 | 2.029 |
| edf10rep2 | 4,115 | 45.8% | 1.542 | 2.294 | 8.540 | 3.699 | 2.232 |
| fig4untreated | 1,010 | 31.4% | 3.103 | 3.298 | 2.541 | 4.102 | 3.255 |
| rtx_rep1 | 2,896 | 23.0% | 3.595 | 4.256 | 6.576 | 3.550 | 3.531 |
| rtx_rep2 | 2,948 | 36.0% | 1.630 | 1.711 | 1.288 | 1.537 | 1.337 |
| rtx_rep3 | 4,125 | 25.9% | 1.466 | 1.670 | 2.559 | 1.717 | 1.478 |
| sdc_corrected | 532 | 2.7% | 3.805 | 4.508 | 4.409 | 3.791 | 3.219 |
| sdc_raw | 379 | 2.0% | 10.029 | 10.945 | 8.363 | 8.889 | 9.189 |

Coverage is relative to the original grouped-code input rows in valid temporal panels, not all localizations in the acquisition. Values have units nm^2; interpreting them as physical variances requires their respective observation models. Quadratic and finite-frequency fits remain different estimands. comparison.json retains raw estimates, both frequencies, Fieller ratio sets, paired spatial-block differences, invalid draws, source rejections, effective sample weights and every historical replay. A confidence interval spanning a difference of zero is not an equivalence result. Spatial blocks are not biological replicates.

## Same-row disagreements and uncertainty

The two-frequency bounded point fit can conceal unsupported ratios. Every frequency-specific difference and paired spatial-block interval is therefore reported separately. Intervals are conditional on the observed oracle selection. Intervals are pointwise and unadjusted for multiple comparisons; no familywise or equivalence claim was specified.

| View | q (nm^-1) | Code minus diagonal equivalent T | Conditional paired pointwise 95% interval | Valid bootstrap draws |
|---|---:|---:|---|---:|
| npc | 0.10 | -0.644 | [-2.268, 0.662] | 100.0% |
| npc | 0.15 | 0.317 | [-0.716, 1.246] | 100.0% |
| edf10rep1 | 0.10 | 0.105 | [-4.670, 2.577] | 100.0% |
| edf10rep1 | 0.15 | 0.103 | [-1.317, 0.898] | 100.0% |
| edf10rep2 | 0.10 | 2.389 | [1.343, 2.970] | 100.0% |
| edf10rep2 | 0.15 | 1.269 | [0.368, 1.882] | 100.0% |
| fig4untreated | 0.10 | 1.418 | [-1.872, 3.942] | 100.0% |
| fig4untreated | 0.15 | 0.722 | [-1.363, 3.085] | 100.0% |
| rtx_rep1 | 0.10 | 0.138 | [-1.650, 2.304] | 100.0% |
| rtx_rep1 | 0.15 | -0.007 | [-0.899, 1.211] | 100.0% |
| rtx_rep2 | 0.10 | 0.316 | [-0.722, 1.150] | 100.0% |
| rtx_rep2 | 0.15 | 0.176 | [-0.464, 1.021] | 100.0% |
| rtx_rep3 | 0.10 | 1.090 | [-0.190, 2.574] | 100.0% |
| rtx_rep3 | 0.15 | 0.064 | [-0.637, 0.733] | 100.0% |
| sdc_corrected | 0.10 | 0.808 | [0.120, 1.510] | 100.0% |
| sdc_corrected | 0.15 | 0.519 | [-0.008, 1.080] | 100.0% |
| sdc_raw | 0.10 | -0.531 | [-1.396, 0.909] | 100.0% |
| sdc_raw | 0.15 | -0.243 | [-0.832, 0.352] | 100.0% |

Unsupported original grouped-code frequency-specific uncertainty:

- edf10rep2, q=0.10: no supported finite equivalent-T interval; Fieller set `unbounded_or_empty`. The bounded two-frequency point fit is not an uncertainty certificate.

## Count, noise and selection conditions

All candidate clusters precede the original source-QC gate. Quantiles below describe the retained empirical ingredients. Within-panel scatter can include merged emitters and motion; it is not independently known localization noise. Full per-panel counts, burst multiplicities, candidate-neighbour distances, frame boundaries and rejection reasons are retained in comparison.json and the NPZ tables.

| View | Candidate clusters | Median candidate → eligible events by panel | Median reported event variance (nm^2) | Median within-panel per-axis scatter (nm^2) | Event rows retained |
|---|---:|---|---:|---:|---:|
| npc | 50,549 | 0/0/0/0 → 4/4/4/4 | 7.092 | 18.982 | 14.0% |
| edf10rep1 | 53,509 | 0/0/0/0 → 4/4/4/4 | 4.400 | 12.023 | 14.5% |
| edf10rep2 | 44,797 | 0/0/0/0 → 6/6/6/6 | 3.741 | 12.108 | 42.7% |
| fig4untreated | 23,591 | 0/0/0/0 → 7/7/7/7 | 5.098 | 26.257 | 28.2% |
| rtx_rep1 | 48,563 | 1/0/0/0 → 5/5/5/5 | 3.820 | 22.797 | 21.6% |
| rtx_rep2 | 32,994 | 0/0/0/0 → 6/6/6/6 | 4.176 | 13.809 | 35.2% |
| rtx_rep3 | 82,777 | 0/0/0/0 → 5/5/5/5 | 10.068 | 26.560 | 24.1% |
| sdc_corrected | 137,093 | 1/0/0/0 → 28.5/29/27/20 | 9.763 | 76.292 | 2.7% |
| sdc_raw | 137,868 | 1/0/0/0 → 29/27/27/21 | 9.763 | 80.217 | 1.9% |

## Preparation-specific reference response

Six CD20 acquisitions complete 96 control worlds and 384,000 group draws: T=0,1,3,5 nm^2; Gaussian and within-panel residual noise; two seeds; 4,000 group draws in every world. Inputs retain all candidate clusters before oracle source QC, their geometry, panel counts, noise and localization multiplicities. Frozen source QC is rerun after adding known window motion. Paired per-group influences preserve common random draws across T within each seed.

| Acquisition | Event-weight slopes (Gaussian / residual) | Event-weight response qualified | Localization-weight slopes (Gaussian / residual) | Localization-weight response qualified |
|---|---:|---|---:|---|
| edf10rep1 | 0.932 / 0.911 | True | 0.897 / 0.872 | True |
| edf10rep2 | 0.926 / 0.981 | True | 0.960 / 0.945 | True |
| fig4untreated | 1.009 / 0.954 | True | 1.045 / 0.946 | False |
| rtx_rep1 | 0.962 / 0.902 | False | 1.010 / 0.847 | False |
| rtx_rep2 | 0.964 / 0.885 | True | 0.902 / 0.915 | False |
| rtx_rep3 | 0.975 / 0.856 | False | 0.969 / 0.841 | False |

Failed criteria, including the unweighted reference sensitivity:

- fig4untreated, localization_count: gaussian maximum affine residual 0.238 >0.2 nm^2; residual maximum affine residual 0.307 >0.2 nm^2.
- rtx_rep1, event_count: gaussian maximum affine residual 0.209 >0.2 nm^2.
- rtx_rep1, localization_count: noise-model slope difference 17.6% >10%; gaussian maximum affine residual 0.261 >0.2 nm^2.
- rtx_rep2, localization_count: gaussian maximum affine residual 0.318 >0.2 nm^2; residual maximum affine residual 0.465 >0.2 nm^2.
- rtx_rep3, event_count: noise-model slope difference 13.0% >10%; gaussian maximum affine residual 0.208 >0.2 nm^2.
- rtx_rep3, localization_count: noise-model slope difference 14.2% >10%; gaussian maximum affine residual 0.428 >0.2 nm^2.

Qualification requires positive slope intervals, <=10% between-model slope difference, and <=0.2 nm^2 maximum affine residual. These are declared sensitivity criteria, not a field standard. Unqualified responses stay in the table and JSON. Intercepts, complete source-selection coverage, inferred-label impurity, model-specific transformed observations and uncertainty are retained. No model was selected to improve code/oracle agreement.

These are field-derived **conditional** controls. Inferred clusters may already merge real emitters; neither noise model establishes physical assignments, event independence or the actual displacement law. Transformed observation intervals do not include every reference-model uncertainty. A qualified affine response cannot certify the acquisition or turn quiet D into validation.

## Temporal support and the existing SDC calibration

NPC panels cover a 12,500-frame acquisition; CD20 panels are quarters of round-specific acquisitions ranging from 20,000 to 50,000 frames; SDC has three 4,200-frame panels and a final 4,250-frame panel, with 50-frame guard gaps over 17,000 frames. These are distinct finite temporal supports, not a universal common lag. Per-source mean event frames, exact localization and event panel counts, and round-specific boundaries are retained. No absolute synchronization of sequential rounds is inferred.

For noise-free panel means a,b,c,d, the quadratic reference equals [(a-c)(b-d)+(a-d)(b-c)]/4. It becomes a per-window variance only under the persistent independent-window model; general motion produces a signed increment-product functional. The finite-q contrast identifies attenuation under its separate assumptions.

The prior SDC simulation failed to reproduce every diagnostic and selection coverage. Its fitted intercept includes modeled settling motion, so subtracting it targets only the injected random component. It cannot be used as a fully matched reference for a code contrast retaining that physical motion. Both raw/corrected observations and the old conditional sensitivity remain. See ESTIMANDS.md.



## Focused reconciliation sensitivity — 23 October 2026


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

The original fitted responses replay for all six CD20 acquisitions. Observation and response-coefficient intervals are propagated by the specified conservative Bonferroni rectangle. Each noise law remains separate; their envelope is a sensitivity envelope, not a confidence interval for every physical model. Calibrated responses apply only to the original QC. NPC and SDC lack these matched response controls and are explicitly unsupported for this transformation.

| Acquisition | Weight | Original qualified | Gaussian propagated interval | Residual propagated interval |
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

## Interpretation of sample size and support

The 81 sensitivity settings reuse nine views of eight acquisitions from two
published studies. They are not 81 independent experiments or a demonstration
of eight independent biological preparations. Gap counts are descriptive.
"Supported" in the multi-source-only diagnostic means eligible finite point-fit
comparison after the declared sign/boundary checks. It does not mean every
individual Fieller set is usable; those topologies are reported separately.
Paired intervals condition on the retained inferred assignments and selection,
without independently propagating physical identity errors. No equivalence,
simultaneous coverage or multiplicity-adjusted discovery is asserted.
