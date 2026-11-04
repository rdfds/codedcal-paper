# Supplementary reporting details

All entries below reaggregate retained reductions without new fits or random draws.

## Full-population arithmetic baselines

All fits have nm²-equivalent units. Historical phases used float32 for CD20 and SDC; NPC already used float64. The centered-float64 full-population baseline is fixed across the nine settings of a view and is exactly the full-code fit already plotted in Figure 2. Selected and diagonal columns use the original eligibility setting (r1-original); neither is physical ground truth. Thus Figure 2 already compares centered-float64 full and selected reductions; the separate earlier historical fits quantify the small arithmetic change, rather than a hidden confound within that display.

| View | Full historical phase convention | Full centered float64 (Figure 2) | Selected float64 | Diagonal float64 |
|---|---:|---:|---:|---:|
| npc | 1.493590 | 1.493590 | 1.051642 | 0.897886 |
| edf10rep1 | 3.681968 | 3.681840 | 2.132884 | 2.029292 |
| edf10rep2 | 8.539803 | 8.539691 | 3.699445 | 2.232026 |
| fig4untreated | 2.540930 | 2.540867 | 4.101551 | 3.255089 |
| rtx_rep1 | 6.576413 | 6.576367 | 3.549755 | 3.531452 |
| rtx_rep2 | 1.288430 | 1.288463 | 1.536903 | 1.337316 |
| rtx_rep3 | 2.558621 | 2.558786 | 1.716770 | 1.478476 |
| sdc_corrected | 4.409094 | 4.409094 | 3.790668 | 3.219131 |
| sdc_raw | 8.363324 | 8.363324 | 8.888647 | 9.189291 |

## Spatial support and validity across every sensitivity setting

The minimum selected-source support is 7 nonempty spatial blocks. The 4-by-3 partition has 12 possible cells; the actual resampling frame uses only tile IDs occupied by the original coarse groups (7–12 cells by view). Cells in that frame with zero selected sources remain in it. The minimum sign-valid draw fraction is 100.00% across all 162 frequency evaluations. A high valid fraction does not imply broad spatial support or calibrated coverage.

The minimum code-input coverage is 0.052831% in sdc_raw, r1.2-stringent, with 11 selected sources and 8 nonempty blocks. It rounds to 0.1% in the retained figure. Event coverage divides selected candidate binding-event rows by all candidate event rows in valid temporal panels; code-input coverage instead divides retained frame-level localization rows by all code-input localization rows in those panels. Burst multiplicity makes these denominators different.

The complete support.csv contains block-level selected-source counts and valid-draw counts for all 81 settings at both frequencies. Three unsupported settings each evaluate two estimators at two frequencies: 3 × 2 × 2 = 12 all-real or disjoint individual Fieller sets.

| View | q | Positive | Negative | Zero-spanning | Unavailable |
|---|---:|---:|---:|---:|---:|
| npc | 0.1 | 0 | 0 | 9 | 0 |
| npc | 0.15 | 5 | 0 | 4 | 0 |
| edf10rep1 | 0.1 | 0 | 0 | 9 | 0 |
| edf10rep1 | 0.15 | 0 | 0 | 9 | 0 |
| edf10rep2 | 0.1 | 4 | 0 | 5 | 0 |
| edf10rep2 | 0.15 | 9 | 0 | 0 | 0 |
| fig4untreated | 0.1 | 1 | 0 | 8 | 0 |
| fig4untreated | 0.15 | 0 | 0 | 9 | 0 |
| rtx_rep1 | 0.1 | 1 | 0 | 8 | 0 |
| rtx_rep1 | 0.15 | 0 | 0 | 9 | 0 |
| rtx_rep2 | 0.1 | 0 | 0 | 9 | 0 |
| rtx_rep2 | 0.15 | 0 | 0 | 9 | 0 |
| rtx_rep3 | 0.1 | 0 | 0 | 9 | 0 |
| rtx_rep3 | 0.15 | 0 | 0 | 9 | 0 |
| sdc_corrected | 0.1 | 1 | 0 | 8 | 0 |
| sdc_corrected | 0.15 | 1 | 0 | 8 | 0 |
| sdc_raw | 0.1 | 0 | 1 | 8 | 0 |
| sdc_raw | 0.15 | 0 | 0 | 9 | 0 |

## Temporal-variant qualification accounting

The inherited four-gate qualification is unchanged. Its usable-and-covering frequency gate is at least 0.90, not 0.95. The final inherited gate asks whether the field-bootstrap range contains nominal 0.95. The companion column only ignores a failure of that last gate from above; it is a descriptive reclassification, not a new qualification criterion or a one-sided 95% confidence guarantee. The union counts overlapping bias, usability and joint-coverage failures once.

| Method | Variant | All inherited gates /144 | Companion allowing above-nominal only /144 | Bias/usability/joint-coverage failure union /144 |
|---|---|---:|---:|---:|
| coded | arithmetic | 115 | 121 | 23 |
| coded | geometric | 2 | 126 | 18 |
| known_source | arithmetic | 99 | 123 | 21 |
| known_source | geometric | 17 | 129 | 15 |
