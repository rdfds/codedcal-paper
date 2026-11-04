# Existing coded-acquisition reference comparison — 7 October 2026

Eight distinct acquisitions / nine raw-or-corrected views are included. Figure4 RTX is the same acquisition as EDF11/EDF12 RTX replicate 1; all four picked-member hashes agree, and it counts once. All analyses are post hoc on exposed public data.

## What is now reconciled

The historical quadratic reference uses weighted binding-event means and event-count products on a clusterable subset. The code statistic uses unweighted localization Fourier sums and localization-count products on the original, broader grouped-code population. A matched assignment-assisted diagonal contrast and code replay on its exact retained localization rows separate these differences. Source selection remains oracle-dependent; this is a validation diagnostic, not a new assignment-free workflow or proof about unselected sources.

All eligible historical source rows, contrasts and count weights are replayed before interpreting new controls. The raw SDC view has no historical raw-oracle source file; its corrected oracle and both full-code moment arrays are audited. Legacy code fits below retain the historical bounded two-frequency fit solely for comparison; they are not unconditional physical variances.

## Interpretation

The reference comparison is now acquisition-complete and target-aware. For RTX replicate 1, the original grouped-code fit is 6.576 nm^2, while same-row code/diagonal fits are 3.550/3.531. Selection accounts for much of that old discrepancy. It does not improve every comparison: Figure4 untreated moves from an original fit of 2.541 to a same-row fit of 4.102, against a diagonal fit of 3.255.

EDF10 replicate 2 retains positive conditional code-minus-diagonal differences at both frequencies; corrected SDC retains one at q=0.10. Other same-row paired intervals span zero, which does not demonstrate equivalence. Oracle selection can affect code exogeneity, and inferred source assignments, heterogeneous motion and changing kernels remain possible limitations. These data cannot uniquely assign the remaining discrepancies to physical motion or estimator bias.

Declared conditional response qualification passes for 4/6 event-weighted and 2/6 localization-weighted CD20 references. An NPC-derived affine response is therefore not a generally validated transfer. The result repairs the validation comparison and its claims; it does not establish physical temporal truth, validation of unselected sources or a new independent benefit over strong alternatives.

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

## Numerical provenance and reproduction

The original CD20 code evaluated phases on global float32 positions, and the original SDC scalar phase operations also retained float32. New target-matching diagnostics use centered float64 positions. reconcile.py re-executes the historical code function to verify its original moments and records the separate numerical sensitivity without loosening replay tolerances. The first three control fields use dense legacy clustering; subsequent controls use a partition-equivalent sparse implementation, checked at the strict cutoff. The initial 400-draw pilot and an import failure before any controls remain; primary controls uniformly use 4,000 draws.

From the repository root, using ../merlin-venv.nosync/bin/python:

```sh
python outputs/reference-reconciliation-oct07/fetch_existing.py
python outputs/reference-reconciliation-oct07/pipeline.py
python outputs/reference-reconciliation-oct07/test_controls.py
python outputs/reference-reconciliation-oct07/reconcile.py
python outputs/reference-reconciliation-oct07/write_report.py
python outputs/reference-reconciliation-oct07/verify.py
```

NPC picked tables and SDC reduced projections are reused from the verified existing caches described in reference_characterize.py. The CD20 fetcher checks whole-member historical SHA-256, two extraction implementations and CRC, retaining one raw copy outside Git. Field-derived simulation ingredients, per-group controls and verification hashes are retained in this package.
