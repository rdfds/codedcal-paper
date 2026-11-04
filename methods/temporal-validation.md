# Supplementary Table 3: Conditional temporal measurement validation

Frozen settings, 9 October 2026. All retained pre-QC candidate clusters from eight public acquisitions / nine views condition synthetic emitters. Source labels are not physically validated. Each noise/view combination has four primary worlds and two pointwise frequencies (eight conditions), 40 simulated fields per world and two code assignments per field. Coarse-group templates are sampled uniformly with replacement; analytic truth averages the full template population. Observed scope preserves group count, not the original composition. The matched functional is a count/kernel-weighted attenuation, not an unweighted molecular variance.

| View | Noise law | Coded gates passed /8 | Known-source gates passed /8 | Coded joint usable coverage range | Median usable interval width range |
|---|---|---:|---:|---:|---:|
| npc | gaussian | 8/8 | 7/8 | 92.5–96.2% | 0.0726–0.0888 |
| npc | residual | 8/8 | 7/8 | 90.0–96.2% | 0.0732–0.0874 |
| edf10rep1 | gaussian | 6/8 | 5/8 | 95.0–98.8% | 0.1266–0.1539 |
| edf10rep1 | residual | 7/8 | 6/8 | 87.5–97.5% | 0.1278–0.1552 |
| edf10rep2 | gaussian | 0/8 | 1/8 | 73.8–92.5% | 0.1194–0.4887 |
| edf10rep2 | residual | 0/8 | 0/8 | 12.5–83.8% | 0.0699–0.2912 |
| fig4untreated | gaussian | 8/8 | 7/8 | 92.5–97.5% | 0.1027–0.1318 |
| fig4untreated | residual | 8/8 | 6/8 | 90.0–97.5% | 0.1047–0.1294 |
| rtx_rep1 | gaussian | 7/8 | 7/8 | 92.5–98.8% | 0.0773–0.1304 |
| rtx_rep1 | residual | 5/8 | 6/8 | 86.2–98.8% | 0.0724–0.1544 |
| rtx_rep2 | gaussian | 7/8 | 4/8 | 92.5–100.0% | 0.0482–0.0966 |
| rtx_rep2 | residual | 8/8 | 6/8 | 91.2–96.2% | 0.0449–0.1230 |
| rtx_rep3 | gaussian | 7/8 | 6/8 | 88.8–97.5% | 0.0676–0.1116 |
| rtx_rep3 | residual | 8/8 | 4/8 | 90.0–96.2% | 0.0650–0.1195 |
| sdc_corrected | gaussian | 6/8 | 7/8 | 87.5–98.8% | 0.0753–0.2093 |
| sdc_corrected | residual | 7/8 | 7/8 | 91.2–100.0% | 0.0250–0.0831 |
| sdc_raw | gaussian | 8/8 | 6/8 | 93.8–96.2% | 0.0777–0.2304 |
| sdc_raw | residual | 7/8 | 7/8 | 88.8–96.2% | 0.0252–0.0802 |

Total coded primary conditions passing all gates: 115/144. Gates: absolute Monte Carlo mean ratio bias ≤0.01; ≥90% bounded positive-ratio confidence calls; ≥90% joint usable physical-target coverage; whole-physical-field percentile bootstrap interval includes nominal 95%. No multiplicity-adjusted certificate is intended. All failures remain in the complete CSV. A bootstrap failure above nominal can be an all-success boundary artifact; conservative post-hoc bounded-field intervals are sensitivities, not substituted gates. More passes do not imply better precision.

Fieller sets rely on a group-mean approximation. A finite interval can still be too broad to resolve motion. Correct coverage and zero-motion exclusion are reported as post-hoc descriptive endpoints. The known-source comparator receives privileged synthetic assignments, not an implementable independent molecular reference. All candidate sources enter generation, but unsupported sources have zero triple weight; effective source/group weight counts and noise-fallback dependence are retained.

Fifty-four stress worlds additionally test window-B kernel changes and static code-dependent geometry. Physical attenuation and actual source-diagonal ratio are separate targets; invalid and quiet D misses remain. These stress results describe the frozen arithmetic estimator, not every existing robustness variant. Sampled historical clustering/QC reports impurity and pure-source kernel-mass coverage; it cannot validate physical source identities.

Source Data 70: `outputs/temporal-validation-oct09/PROTOCOL.md`, `freeze.json`, `summary.csv`, `summary.json`, `population-information.json`, `verification.json`, `artifact-manifest.json`, and all per-world truth/group observations. All 126 worlds / 3,528 simulated field realizations complete. Source Data 69 supplies the hashed exposed event templates. Raw/corrected SDC are paired views; Figure4 RTX aliases RTX replicate 1 and is counted once.
