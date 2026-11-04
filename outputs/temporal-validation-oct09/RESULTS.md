# Preparation-derived temporal validation

Completed 10 October 2026; settings frozen 9 October. Complete, frozen **conditional simulation** from eight public acquisitions / nine raw-corrected views. No new acquisition, physical randomized-readout intervention or independent physical source truth.

All 126 prescribed worlds completed: 72 primary and 54 stress worlds, 3,528 simulated physical-field realizations, 7,056 coded observations, and 7,662,816 group/code draws. Paired codes share a physical field; they are clustered in the bootstrap. All 42,336 pointwise cases, 126 analytic truth tables, nine source catalogs and retained artifact hashes passed recount. This audit pass is not a scientific qualification pass.

## Main finding

At the observed group counts, **115/144 coded primary world/frequency conditions pass all frozen gates**, compared with **99/144** for the source-label-informed Fourier comparator; **89** pass for both. Gates require absolute mean ratio bias ≤0.01, ≥90% finite positive confidence calls, ≥90% joint usable physical-target coverage, and a whole-field bootstrap interval including nominal 95%. These are pointwise study criteria, not established field standards or multiplicity-adjusted certification.

Fieller sets use a group-mean normal approximation with a Student critical value; they are nominal confidence sets, not exact finite-sample coverage guarantees. Highly concentrated count/kernel weights can leave far fewer effectively contributing groups than the raw group count. Coverage failures in either method remain failures. More gate passes do not imply better precision: broader intervals may cover more often.

The group-common world shares motion within one synthetic pore/group and draws it independently across groups. A motion component shared by an entire acquisition is outside this independent-group uncertainty test. The persistent-window law is imposed, not fitted from real frame lags.

The source-informed comparator has privileged synthetic labels and is an information benchmark. It does not demonstrate an implementable stronger pipeline, and this experiment does not establish a benefit over Shaw–Veatch, which measures a different reference-relative lag quantity. No optimizer boundary or quiet diagnostic certifies a measurement.

Coded primary gate failures (overlap allowed): {'nominal_coverage_in_bootstrap': 20, 'usable_coverage': 21, 'ratio_bias': 14, 'usable_fraction': 13}. Failed criteria remain failures; uncertainty intervals and extreme ratio tails are retained rather than discarded.

Nominal bootstrap failures by direction: {'none': 124, 'above_nominal': 7, 'below_nominal': 13}. A percentile bootstrap at all-success coverage is degenerate; above-nominal failures are not evidence of undercoverage. A separately labeled post-hoc Chernoff interval for the mean bounded physical-field coverage is included for every condition, preserving arbitrary dependence between its two code assignments. It is a finite-sample sensitivity and does not replace any frozen gate.

EDF10 replicate 2 is the strongest failure case: the coded finite-ratio sample mean has maximum absolute bias 0.407, minimum joint usable coverage 12.5%, and a ratio tail reaching 23.7. The known-source comparator reaches only 5% joint coverage in its worst heterogeneous residual condition. Median agreement or a correct extraction identity cannot erase these failures. Finite-ratio mean/RMSE diagnostics are sample summaries; finite population ratio moments are not assumed.

## All views and noise laws

| View | Noise | Coded qualified /8 | Known-source /8 | Coded usable range | Coded joint coverage range | Worst absolute ratio bias |
|---|---|---:|---:|---:|---:|---:|
| npc | gaussian | 8/8 | 7/8 | 100.0%–100.0% | 92.5%–96.2% | 0.0038 |
| npc | residual | 8/8 | 7/8 | 100.0%–100.0% | 90.0%–96.2% | 0.0020 |
| edf10rep1 | gaussian | 6/8 | 5/8 | 100.0%–100.0% | 95.0%–98.8% | 0.0071 |
| edf10rep1 | residual | 7/8 | 6/8 | 100.0%–100.0% | 87.5%–97.5% | 0.0069 |
| edf10rep2 | gaussian | 0/8 | 1/8 | 82.5%–92.5% | 73.8%–92.5% | 0.4071 |
| edf10rep2 | residual | 0/8 | 0/8 | 51.2%–86.2% | 12.5%–83.8% | 0.1300 |
| fig4untreated | gaussian | 8/8 | 7/8 | 100.0%–100.0% | 92.5%–97.5% | 0.0080 |
| fig4untreated | residual | 8/8 | 6/8 | 100.0%–100.0% | 90.0%–97.5% | 0.0085 |
| rtx_rep1 | gaussian | 7/8 | 7/8 | 100.0%–100.0% | 92.5%–98.8% | 0.0113 |
| rtx_rep1 | residual | 5/8 | 6/8 | 100.0%–100.0% | 86.2%–98.8% | 0.0059 |
| rtx_rep2 | gaussian | 7/8 | 4/8 | 100.0%–100.0% | 92.5%–100.0% | 0.0079 |
| rtx_rep2 | residual | 8/8 | 6/8 | 100.0%–100.0% | 91.2%–96.2% | 0.0042 |
| rtx_rep3 | gaussian | 7/8 | 6/8 | 100.0%–100.0% | 88.8%–97.5% | 0.0053 |
| rtx_rep3 | residual | 8/8 | 4/8 | 100.0%–100.0% | 90.0%–96.2% | 0.0087 |
| sdc_corrected | gaussian | 6/8 | 7/8 | 100.0%–100.0% | 87.5%–98.8% | 0.0095 |
| sdc_corrected | residual | 7/8 | 7/8 | 100.0%–100.0% | 91.2%–100.0% | 0.0023 |
| sdc_raw | gaussian | 8/8 | 6/8 | 100.0%–100.0% | 93.8%–96.2% | 0.0080 |
| sdc_raw | residual | 7/8 | 7/8 | 100.0%–100.0% | 88.8%–96.2% | 0.0048 |

These are repeated simulations conditioned on exposed preparations, not eighteen independent acquisitions. SDC raw/corrected share one acquisition; the Figure4 RTX alias is counted once. The 512-group sensitivity uses the same draws and is not independent confirmation.

## Measurement resolution

Finite positive Fieller uncertainty does not imply useful precision for small motion. The table reports pointwise exclusion of zero motion (attenuation ratio=1), plus correct target coverage; these descriptive endpoints were added after early outcomes and do not change frozen qualification.

| View | Noise | Homogeneous T=3: correct coverage and zero exclusion q=.10 / .15 | Median interval width q=.10 / .15 | Known-source correct zero exclusion q=.10 / .15 |
|---|---|---:|---:|---:|
| npc | gaussian | 25.0% / 78.8% | 0.0742 / 0.0855 | 97.5% / 92.5% |
| npc | residual | 31.2% / 85.0% | 0.0737 / 0.0860 | 97.5% / 97.5% |
| edf10rep1 | gaussian | 10.0% / 41.2% | 0.1296 / 0.1539 | 97.5% / 97.5% |
| edf10rep1 | residual | 21.2% / 32.5% | 0.1286 / 0.1521 | 95.0% / 97.5% |
| edf10rep2 | gaussian | 22.5% / 33.8% | 0.1238 / 0.1385 | 45.0% / 45.0% |
| edf10rep2 | residual | 18.8% / 32.5% | 0.0979 / 0.1594 | 30.0% / 32.5% |
| fig4untreated | gaussian | 10.0% / 52.5% | 0.1048 / 0.1278 | 97.5% / 92.5% |
| fig4untreated | residual | 22.5% / 51.2% | 0.1056 / 0.1287 | 97.5% / 95.0% |
| rtx_rep1 | gaussian | 40.0% / 50.0% | 0.0773 / 0.1182 | 97.5% / 100.0% |
| rtx_rep1 | residual | 31.2% / 46.2% | 0.0794 / 0.1165 | 92.5% / 90.0% |
| rtx_rep2 | gaussian | 60.0% / 88.8% | 0.0516 / 0.0784 | 97.5% / 100.0% |
| rtx_rep2 | residual | 63.7% / 86.2% | 0.0473 / 0.0701 | 97.5% / 97.5% |
| rtx_rep3 | gaussian | 30.0% / 52.5% | 0.0699 / 0.1075 | 100.0% / 100.0% |
| rtx_rep3 | residual | 32.5% / 65.0% | 0.0652 / 0.0956 | 100.0% / 100.0% |
| sdc_corrected | gaussian | 35.0% / 21.2% | 0.0753 / 0.1999 | 95.0% / 75.0% |
| sdc_corrected | residual | 93.8% / 92.5% | 0.0272 / 0.0419 | 97.5% / 95.0% |
| sdc_raw | gaussian | 27.5% / 17.5% | 0.0785 / 0.2225 | 95.0% / 67.5% |
| sdc_raw | residual | 95.0% / 88.8% | 0.0272 / 0.0423 | 97.5% / 92.5% |

## Physical interpretation and diagnostic misses

Violations test the frozen arithmetic ratio. The existing split geometric construction can cancel source-uniform multiplicative kernel changes under its own assumptions; these results must not be described as failures of every temporal variant. No new strongest-variant comparison was frozen here.

The physical attenuation target and observable source-diagonal ratio are recorded separately in every world. Window-B noise changes can preserve randomized source extraction while changing the ratio’s interpretation. Static code-dependent positions can break extraction itself. Heterogeneous source motion targets the exact count/kernel-weighted attenuation, not arithmetic motion variance.

| Violation | Coded cases | Errors >0.01 | Such errors with quiet D | Invalid D | Alerts / valid D | Observable-to-physical target difference range |
|---|---:|---:|---:|---:|---:|---:|
| noise_shift | 864 | 727 | 57 | 2 | 92.3% | 0.0043–0.0786 |
| source_noise_shift | 864 | 841 | 26 | 6 | 97.0% | 0.0236–0.6000 |
| code_offset | 864 | 539 | 480 | 2 | 8.2% | 0.0000–0.0000 |

Across common-motion null coded cases, D produces 509/8630 valid-test alerts and 10 invalid tests. These pooled counts are descriptive; acquisitions, frequencies and paired codes are not independent biological replicates. Heterogeneous source motion can produce different weighted window targets despite stable kernels; its alerts are not false alarms under a common-motion null. Both analytic window D targets and per-condition rates are in `summary.json`.

## Sampled reference coverage

Only homogeneous T=3 uses the frozen first-eight-fields, first-sixteen-groups clustering sample. Historical radii/QC remain fixed. Mixed assignments earn no known-source coverage; fragments receive their captured count-triple fraction. This diagnostic tests inferred synthetic identities, not independent biological labeling.

| View | Noise | Passing clusters | Passing impurity | Localization coverage | Pure stationary kernel mass retained |
|---|---|---:|---:|---:|---:|
| npc | gaussian | 280 | 23.6% | 13.8% | 19.7% |
| npc | residual | 278 | 24.1% | 13.0% | 20.1% |
| edf10rep1 | gaussian | 413 | 31.0% | 15.9% | 29.9% |
| edf10rep1 | residual | 331 | 35.0% | 13.3% | 26.3% |
| edf10rep2 | gaussian | 997 | 22.1% | 47.5% | 48.8% |
| edf10rep2 | residual | 889 | 19.0% | 42.3% | 30.2% |
| fig4untreated | gaussian | 460 | 35.7% | 29.5% | 27.3% |
| fig4untreated | residual | 479 | 34.2% | 32.5% | 32.9% |
| rtx_rep1 | gaussian | 575 | 42.3% | 22.6% | 13.4% |
| rtx_rep1 | residual | 557 | 45.1% | 21.0% | 14.2% |
| rtx_rep2 | gaussian | 720 | 31.9% | 36.3% | 15.4% |
| rtx_rep2 | residual | 621 | 30.4% | 34.2% | 20.0% |
| rtx_rep3 | gaussian | 543 | 33.0% | 24.9% | 24.8% |
| rtx_rep3 | residual | 562 | 42.3% | 25.7% | 20.7% |
| sdc_corrected | gaussian | 99 | 70.7% | 1.5% | 0.3% |
| sdc_corrected | residual | 170 | 60.6% | 3.4% | 1.6% |
| sdc_raw | gaussian | 93 | 55.9% | 1.6% | 1.3% |
| sdc_raw | residual | 135 | 63.7% | 2.4% | 0.3% |

## What this changes

This closes the missing preparation-derived finite-observation test under its declared laws. It identifies exactly which conditions support useful uncertainty and which do not, while retaining all source/QC failures. It does **not** establish that the acquired source assignments or persistent-window law describe the instrument, validate an unperformed randomized MERFISH intervention, or resolve the independent molecular-benefit gap. Passing conditions support a conditional measurement claim; failed conditions limit it.

No favorable subset becomes the headline for unrestricted accuracy. The complete table contains every preparation, noise model, frequency, scope and comparator; the summary plot shows all primary observed-count conditions. All original failed S3/ISS/molecular results and density qualifications remain unchanged.

## Retained evidence

`PROTOCOL.md` and `freeze.json` predate outcomes. `controls.json` records analytic extraction, kernel, persistent-motion and Fieller checks; `IMPLEMENTATION.md` records the pre-outcome type-label repair. `catalog-*.npz` holds all retained inferred sources; hashed empirical inputs hold event residuals. Each `worlds/*/truth.npz` holds physical/observable expected moments; every `field-*.npz` retains group observations and template IDs, with seeds/uncertainty/diagnostics/QC in its JSON. `summary.csv` is the full condition table; `verification.json` and `artifact-manifest.json` provide recount and hashes. No outcome-driven additional simulation or threshold selection was performed.


## What population is measured

All retained candidate clusters enter generation, but sources lacking three-panel support have zero diagonal weight. The observable is count/kernel weighted, not the mean budget of every candidate or raw localization. Noise fallback is inferred from the preparation; it is not measured physical truth. Effective weight counts below quantify concentration, not independent biological replication.

| View | Noise | Candidates | Triple support | Effective sources q=.10 / .15 | Effective groups q=.10 / .15 | Heterogeneous apparent T q=.10 / .15 (source mean=3) |
|---|---|---:|---:|---:|---:|---:|
| npc | gaussian | 50549 | 14203 | 813.9 / 911.2 | 498.6 / 529.1 | 4.99 / 4.20 |
| npc | residual | 50549 | 14203 | 802.6 / 873.8 | 494.0 / 515.9 | 4.99 / 4.36 |
| edf10rep1 | gaussian | 53509 | 17034 | 1453.2 / 1733.9 | 643.5 / 690.1 | 5.19 / 4.44 |
| edf10rep1 | residual | 53509 | 17034 | 1469.7 / 1748.3 | 646.5 / 692.4 | 5.14 / 4.45 |
| edf10rep2 | gaussian | 44797 | 10457 | 536.4 / 1191.2 | 244.7 / 540.6 | 8.31 / 4.98 |
| edf10rep2 | residual | 44797 | 10457 | 6.5 / 11.0 | 2.9 / 4.5 | 93.53 / 34.58 |
| fig4untreated | gaussian | 23591 | 2918 | 327.7 / 299.7 | 197.3 / 186.6 | 5.38 / 4.28 |
| fig4untreated | residual | 23591 | 2918 | 323.9 / 298.3 | 196.5 / 186.9 | 5.55 / 4.75 |
| rtx_rep1 | gaussian | 48563 | 11868 | 584.8 / 798.5 | 341.4 / 435.1 | 10.69 / 6.21 |
| rtx_rep1 | residual | 48563 | 11868 | 478.9 / 471.6 | 278.5 / 240.5 | 12.83 / 9.97 |
| rtx_rep2 | gaussian | 32994 | 8926 | 600.8 / 850.7 | 343.4 / 438.4 | 12.18 / 7.62 |
| rtx_rep2 | residual | 32994 | 8926 | 486.2 / 606.4 | 268.2 / 329.3 | 15.00 / 11.12 |
| rtx_rep3 | gaussian | 82777 | 18665 | 1435.1 / 1694.7 | 567.8 / 671.3 | 7.45 / 5.16 |
| rtx_rep3 | residual | 82777 | 18665 | 1153.2 / 1254.3 | 471.7 / 502.9 | 9.18 / 7.49 |
| sdc_corrected | gaussian | 137093 | 13034 | 634.0 / 392.2 | 412.2 / 290.2 | 9.95 / 6.39 |
| sdc_corrected | residual | 137093 | 13034 | 766.9 / 701.2 | 467.1 / 437.4 | 21.77 / 20.23 |
| sdc_raw | gaussian | 137868 | 12809 | 620.2 / 383.4 | 411.4 / 290.4 | 9.62 / 6.08 |
| sdc_raw | residual | 137868 | 12809 | 768.2 / 701.8 | 470.0 / 439.0 | 21.07 / 19.68 |