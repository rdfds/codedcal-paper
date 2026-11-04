# Supplementary Table 4: existing temporal variants

Source Data 71 reuses all 126 worlds and 3,528 simulated fields; 84,672 cases
are independently recounted. This is conditional computation, not physical
validation. Original Source Data 70 results and failed gates remain unchanged.



The existing geometric construction cancels a source-uniform multiplicative
kernel change under the declared common-motion law. The arithmetic construction
does not generally cancel it. On the specified uniform-noise stress, median
absolute Monte Carlo mean bias across 36 coded view/noise/frequency conditions
falls from 0.0292 to 0.0043; mean joint usable coverage rises from 72.6% to 96.9%.
This is a conditional robustness benefit, not physical instrument validation.

Its conservative joint sets are wider and often lose small-motion resolution.
For NPC Gaussian T=3, correct coverage with positive-motion resolution falls
from 25.0% to 0% at q=.10 and from 78.8% to 11.2% at q=.15. EDF10 replicate 2
retains poor support. Source-dependent kernel and code-geometry violations remain;
quiet-D failures are not repaired. No data-selected switch or safe automatic
recommendation is validated.

Arithmetic uses a pointwise 95% Fieller set. Geometric uses two 97.5% window
Fieller sets propagated through their positive product, yielding a conservative
nominal 95% joint rectangle when marginal approximations are calibrated. No
independence between windows is assumed; both require independent groups.
Concentrated contributions can defeat either approximation. More coverage alone
is not more information. Gaussian-equivalent variances are conditional on the
motion model, not total molecular-distance errors.

## All declared violations

Each stress has 864 observed-count coded cases. Coverage averages and median
absolute *sample-mean bias* are across conditions, not independent preparations.
Available ratios and errors do not imply usable confidence support.

| Stress | Variant | Available cases | Mean joint usable coverage | Median absolute sample-mean ratio bias | Errors >0.01 / quiet-D errors |
|---|---|---:|---:|---:|---:|
| noise_shift | arithmetic | 863/864 | 72.6% | 0.0292 | 727 / 57 |
| noise_shift | geometric | 862/864 | 96.9% | 0.0043 | 601 / 57 |
| source_noise_shift | arithmetic | 864/864 | 29.9% | 0.1172 | 841 / 26 |
| source_noise_shift | geometric | 858/864 | 83.4% | 0.0415 | 772 / 26 |
| code_offset | arithmetic | 864/864 | 92.6% | 0.0039 | 539 / 480 |
| code_offset | geometric | 862/864 | 94.9% | 0.0036 | 543 / 490 |

## Availability-matched diagnostic (post hoc)

Point availability differs between variants. Recounting identical available
paired draws retains 862/864 uniform-kernel stress cases and the same median
absolute sample-mean biases, 0.0292 versus 0.0043. Source-dependent stress
retains 858/864, with 0.1267 versus 0.0415; code geometry retains 862/864,
with 0.0038 versus 0.0036. This post-hoc fairness diagnostic is not a replacement
gate or an estimate over excluded draws. Unconditional coverage denominators
remain all 864 cases. Full records are in availability-matched.json.

## Inherited gate accounting, not a ranking of conservative intervals

The original four gates are retained literally. The nominal-in-bootstrap gate
asks whether coverage is compatible with exactly 95%. That is not an appropriate
accuracy ranking for deliberately conservative joint rectangles: 100% observed
coverage can fail it from above. Neither those failures nor all-success bootstrap
degeneracy demonstrate undercoverage. Do not advertise geometric pass counts as
evidence that it is universally worse, or its greater coverage as a universal win.
Every failure, including below-nominal coverage and bias, remains below and in CSV.

| Method | Variant | All inherited gates /144 | Above / below nominal bootstrap failures | Bias / usability / joint-coverage failures |
|---|---|---:|---:|---:|
| coded | arithmetic | 115/144 | 7 / 13 | 14 / 13 / 21 |
| coded | geometric | 2/144 | 126 / 16 | 15 / 16 / 16 |
| known_source | arithmetic | 99/144 | 25 / 7 | 10 / 6 / 15 |
| known_source | geometric | 17/144 | 112 / 10 | 11 / 11 / 11 |

## Targets and independent recount

Each variant is scored against its own baseline physical functional, and its
actual-diagonal observable target is recorded separately. Under heterogeneous
source motion the geometric and arithmetic targets can differ: changing the
weighted functional is not recovering the original arithmetic target. summary.csv
retains the original-target error and target change as well. The known-source
comparator has privileged synthetic labels and is not a deployable alternative.

All 42,336 arithmetic cases and 1,008 summaries reproduce Source Data 70,
including its 115/144 coded qualification verdict. Independent frozen-formula
recounts check both geometric windows in 84,672 comparisons. The initial
verification stop is retained: raw versus normalized float64 summation differs
by 1.44e-12 relatively for an already-unusable, nearly cancelled reference.
verify_v2 uses the pre-existing arithmetic point-replay tolerance of 1e-10;
observations, confidence sets, summary numbers and scientific gates did not change.
The original verifier, failed log and precision-only amendment remain.

Reproduction (repository root, pinned analysis Python): controls.py, run.py freeze
in a new working copy before outcomes; run.py execute; run.py summarize;
verify_v2.py; report.py. The existing freeze and completed outcomes refuse silent
replacement. All source hashes and compressed per-case sets are retained.

## Complete observed-count coded primary table

| View | Noise | q | World | Variant | Bias | Usable | Joint coverage | Width | Correct coverage + positive-motion resolution | Gates |
|---|---|---:|---|---|---:|---:|---:|---:|---:|---|
| edf10rep1 | gaussian | 0.10 | group_common | arithmetic | -0.0007 | 100.0% | 96.2% | 0.1285 | 12.5% | all pass |
| edf10rep1 | gaussian | 0.10 | group_common | geometric | -0.0005 | 100.0% | 100.0% | 0.2552 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | gaussian | 0.15 | group_common | arithmetic | 0.0020 | 100.0% | 95.0% | 0.1483 | 32.5% | all pass |
| edf10rep1 | gaussian | 0.15 | group_common | geometric | 0.0023 | 100.0% | 100.0% | 0.2897 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | gaussian | 0.10 | heterogeneous | arithmetic | -0.0027 | 100.0% | 97.5% | 0.1289 | 32.5% | all pass |
| edf10rep1 | gaussian | 0.10 | heterogeneous | geometric | -0.0025 | 100.0% | 100.0% | 0.2548 | 1.2% | nominal_coverage_in_bootstrap |
| edf10rep1 | gaussian | 0.15 | heterogeneous | arithmetic | -0.0005 | 100.0% | 98.8% | 0.1495 | 68.8% | nominal_coverage_in_bootstrap |
| edf10rep1 | gaussian | 0.15 | heterogeneous | geometric | -0.0001 | 100.0% | 100.0% | 0.2836 | 7.5% | nominal_coverage_in_bootstrap |
| edf10rep1 | gaussian | 0.10 | homogeneous | arithmetic | 0.0015 | 100.0% | 97.5% | 0.1296 | 10.0% | all pass |
| edf10rep1 | gaussian | 0.10 | homogeneous | geometric | 0.0021 | 100.0% | 100.0% | 0.2565 | 1.2% | nominal_coverage_in_bootstrap |
| edf10rep1 | gaussian | 0.15 | homogeneous | arithmetic | 0.0020 | 100.0% | 97.5% | 0.1539 | 41.2% | all pass |
| edf10rep1 | gaussian | 0.15 | homogeneous | geometric | 0.0027 | 100.0% | 100.0% | 0.2928 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | gaussian | 0.10 | zero | arithmetic | -0.0071 | 100.0% | 97.5% | 0.1266 | 0.0% | all pass |
| edf10rep1 | gaussian | 0.10 | zero | geometric | -0.0068 | 100.0% | 100.0% | 0.2579 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | gaussian | 0.15 | zero | arithmetic | -0.0061 | 100.0% | 98.8% | 0.1529 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | gaussian | 0.15 | zero | geometric | -0.0054 | 100.0% | 100.0% | 0.3007 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | residual | 0.10 | group_common | arithmetic | 0.0018 | 100.0% | 95.0% | 0.1285 | 7.5% | all pass |
| edf10rep1 | residual | 0.10 | group_common | geometric | 0.0025 | 100.0% | 100.0% | 0.2655 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | residual | 0.15 | group_common | arithmetic | 0.0018 | 100.0% | 95.0% | 0.1505 | 32.5% | all pass |
| edf10rep1 | residual | 0.15 | group_common | geometric | 0.0025 | 100.0% | 100.0% | 0.2897 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | residual | 0.10 | heterogeneous | arithmetic | 0.0069 | 100.0% | 87.5% | 0.1289 | 23.8% | usable_coverage |
| edf10rep1 | residual | 0.10 | heterogeneous | geometric | 0.0073 | 100.0% | 100.0% | 0.2589 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | residual | 0.15 | heterogeneous | arithmetic | 0.0041 | 100.0% | 92.5% | 0.1514 | 53.8% | all pass |
| edf10rep1 | residual | 0.15 | heterogeneous | geometric | 0.0047 | 100.0% | 100.0% | 0.2888 | 6.2% | nominal_coverage_in_bootstrap |
| edf10rep1 | residual | 0.10 | homogeneous | arithmetic | 0.0024 | 100.0% | 97.5% | 0.1286 | 21.2% | all pass |
| edf10rep1 | residual | 0.10 | homogeneous | geometric | 0.0028 | 100.0% | 100.0% | 0.2598 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | residual | 0.15 | homogeneous | arithmetic | 0.0033 | 100.0% | 93.8% | 0.1521 | 32.5% | all pass |
| edf10rep1 | residual | 0.15 | homogeneous | geometric | 0.0039 | 100.0% | 100.0% | 0.2915 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | residual | 0.10 | zero | arithmetic | 0.0029 | 100.0% | 96.2% | 0.1278 | 0.0% | all pass |
| edf10rep1 | residual | 0.10 | zero | geometric | 0.0034 | 100.0% | 100.0% | 0.2603 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep1 | residual | 0.15 | zero | arithmetic | 0.0040 | 100.0% | 92.5% | 0.1552 | 0.0% | all pass |
| edf10rep1 | residual | 0.15 | zero | geometric | 0.0042 | 100.0% | 100.0% | 0.3022 | 0.0% | nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.10 | group_common | arithmetic | -0.0170 | 92.5% | 92.5% | 0.2707 | 15.0% | ratio_bias |
| edf10rep2 | gaussian | 0.10 | group_common | geometric | 0.0066 | 71.2% | 71.2% | 0.4589 | 0.0% | usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.15 | group_common | arithmetic | 0.0529 | 88.8% | 88.8% | 0.4887 | 32.5% | ratio_bias,usable_fraction,usable_coverage |
| edf10rep2 | gaussian | 0.15 | group_common | geometric | 0.6620 | 71.2% | 71.2% | 0.5893 | 3.8% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.10 | heterogeneous | arithmetic | 0.1293 | 82.5% | 73.8% | 0.1443 | 37.5% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.10 | heterogeneous | geometric | 0.1499 | 71.2% | 71.2% | 0.2179 | 16.2% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.15 | heterogeneous | arithmetic | 0.1022 | 85.0% | 80.0% | 0.1194 | 41.2% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.15 | heterogeneous | geometric | 0.1851 | 67.5% | 67.5% | 0.1615 | 35.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.10 | homogeneous | arithmetic | 0.0135 | 90.0% | 85.0% | 0.1238 | 22.5% | ratio_bias,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.10 | homogeneous | geometric | 0.0498 | 76.2% | 76.2% | 0.1914 | 1.2% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.15 | homogeneous | arithmetic | 0.0869 | 83.8% | 82.5% | 0.1385 | 33.8% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.15 | homogeneous | geometric | 0.0204 | 67.5% | 67.5% | 0.1676 | 15.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.10 | zero | arithmetic | 0.0539 | 90.0% | 88.8% | 0.3159 | 0.0% | ratio_bias,usable_coverage |
| edf10rep2 | gaussian | 0.10 | zero | geometric | 0.0727 | 77.5% | 77.5% | 0.4579 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.15 | zero | arithmetic | 0.4071 | 83.8% | 82.5% | 0.3128 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | gaussian | 0.15 | zero | geometric | 0.1896 | 66.2% | 66.2% | 0.4695 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.10 | group_common | arithmetic | 0.0024 | 66.2% | 66.2% | 0.1117 | 17.5% | usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.10 | group_common | geometric | 0.0147 | 56.2% | 56.2% | 0.2382 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.15 | group_common | arithmetic | -0.0180 | 80.0% | 76.2% | 0.1797 | 31.2% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.15 | group_common | geometric | 0.0020 | 55.0% | 55.0% | 0.2045 | 5.0% | usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.10 | heterogeneous | arithmetic | 0.1300 | 51.2% | 12.5% | 0.2144 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.10 | heterogeneous | geometric | 0.2402 | 32.5% | 3.8% | 0.1931 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.15 | heterogeneous | arithmetic | 0.0774 | 61.3% | 22.5% | 0.2912 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.15 | heterogeneous | geometric | 0.1380 | 43.8% | 8.8% | 0.2562 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.10 | homogeneous | arithmetic | 0.0070 | 62.5% | 62.5% | 0.0979 | 18.8% | usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.10 | homogeneous | geometric | 0.0239 | 50.0% | 50.0% | 0.1615 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.15 | homogeneous | arithmetic | 0.0041 | 78.8% | 78.8% | 0.1594 | 32.5% | usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.15 | homogeneous | geometric | 0.0368 | 53.8% | 53.8% | 0.1700 | 7.5% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.10 | zero | arithmetic | 0.0312 | 71.2% | 67.5% | 0.0699 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.10 | zero | geometric | 0.0422 | 61.3% | 61.3% | 0.1616 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.15 | zero | arithmetic | -0.0222 | 86.2% | 83.8% | 0.1193 | 0.0% | ratio_bias,usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| edf10rep2 | residual | 0.15 | zero | geometric | 0.0060 | 62.5% | 62.5% | 0.1727 | 0.0% | usable_fraction,usable_coverage,nominal_coverage_in_bootstrap |
| fig4untreated | gaussian | 0.10 | group_common | arithmetic | -0.0049 | 100.0% | 96.2% | 0.1027 | 13.8% | all pass |
| fig4untreated | gaussian | 0.10 | group_common | geometric | -0.0047 | 100.0% | 100.0% | 0.2029 | 1.2% | nominal_coverage_in_bootstrap |
| fig4untreated | gaussian | 0.15 | group_common | arithmetic | -0.0080 | 100.0% | 96.2% | 0.1292 | 55.0% | all pass |
| fig4untreated | gaussian | 0.15 | group_common | geometric | -0.0075 | 100.0% | 100.0% | 0.2455 | 3.8% | nominal_coverage_in_bootstrap |
| fig4untreated | gaussian | 0.10 | heterogeneous | arithmetic | -0.0014 | 100.0% | 96.2% | 0.1055 | 36.2% | all pass |
| fig4untreated | gaussian | 0.10 | heterogeneous | geometric | -0.0009 | 100.0% | 100.0% | 0.2080 | 3.8% | nominal_coverage_in_bootstrap |
| fig4untreated | gaussian | 0.15 | heterogeneous | arithmetic | 0.0001 | 100.0% | 95.0% | 0.1311 | 71.2% | all pass |
| fig4untreated | gaussian | 0.15 | heterogeneous | geometric | 0.0007 | 100.0% | 100.0% | 0.2486 | 6.2% | nominal_coverage_in_bootstrap |
| fig4untreated | gaussian | 0.10 | homogeneous | arithmetic | -0.0004 | 100.0% | 97.5% | 0.1048 | 10.0% | all pass |
| fig4untreated | gaussian | 0.10 | homogeneous | geometric | 0.0001 | 100.0% | 100.0% | 0.2049 | 0.0% | nominal_coverage_in_bootstrap |
| fig4untreated | gaussian | 0.15 | homogeneous | arithmetic | 0.0009 | 100.0% | 96.2% | 0.1278 | 52.5% | all pass |
| fig4untreated | gaussian | 0.15 | homogeneous | geometric | 0.0014 | 100.0% | 100.0% | 0.2526 | 1.2% | nominal_coverage_in_bootstrap |
| fig4untreated | gaussian | 0.10 | zero | arithmetic | 0.0005 | 100.0% | 97.5% | 0.1036 | 0.0% | all pass |
| fig4untreated | gaussian | 0.10 | zero | geometric | 0.0010 | 100.0% | 100.0% | 0.2066 | 0.0% | nominal_coverage_in_bootstrap |
| fig4untreated | gaussian | 0.15 | zero | arithmetic | 0.0006 | 100.0% | 92.5% | 0.1318 | 0.0% | all pass |
| fig4untreated | gaussian | 0.15 | zero | geometric | 0.0016 | 100.0% | 100.0% | 0.2628 | 0.0% | nominal_coverage_in_bootstrap |
| fig4untreated | residual | 0.10 | group_common | arithmetic | -0.0029 | 100.0% | 96.2% | 0.1048 | 17.5% | all pass |
| fig4untreated | residual | 0.10 | group_common | geometric | -0.0023 | 100.0% | 100.0% | 0.1997 | 0.0% | nominal_coverage_in_bootstrap |
| fig4untreated | residual | 0.15 | group_common | arithmetic | -0.0060 | 100.0% | 93.8% | 0.1232 | 55.0% | all pass |
| fig4untreated | residual | 0.15 | group_common | geometric | -0.0051 | 100.0% | 100.0% | 0.2399 | 6.2% | nominal_coverage_in_bootstrap |
| fig4untreated | residual | 0.10 | heterogeneous | arithmetic | 0.0058 | 100.0% | 90.0% | 0.1047 | 47.5% | all pass |
| fig4untreated | residual | 0.10 | heterogeneous | geometric | 0.0069 | 100.0% | 100.0% | 0.1993 | 2.5% | nominal_coverage_in_bootstrap |
| fig4untreated | residual | 0.15 | heterogeneous | arithmetic | 0.0074 | 100.0% | 92.5% | 0.1294 | 76.2% | all pass |
| fig4untreated | residual | 0.15 | heterogeneous | geometric | 0.0078 | 100.0% | 100.0% | 0.2390 | 20.0% | nominal_coverage_in_bootstrap |
| fig4untreated | residual | 0.10 | homogeneous | arithmetic | -0.0003 | 100.0% | 96.2% | 0.1056 | 22.5% | all pass |
| fig4untreated | residual | 0.10 | homogeneous | geometric | -0.0002 | 100.0% | 100.0% | 0.2037 | 0.0% | nominal_coverage_in_bootstrap |
| fig4untreated | residual | 0.15 | homogeneous | arithmetic | 0.0024 | 100.0% | 97.5% | 0.1287 | 51.2% | all pass |
| fig4untreated | residual | 0.15 | homogeneous | geometric | 0.0031 | 100.0% | 100.0% | 0.2425 | 1.2% | nominal_coverage_in_bootstrap |
| fig4untreated | residual | 0.10 | zero | arithmetic | 0.0034 | 100.0% | 93.8% | 0.1071 | 0.0% | all pass |
| fig4untreated | residual | 0.10 | zero | geometric | 0.0039 | 100.0% | 100.0% | 0.2075 | 0.0% | nominal_coverage_in_bootstrap |
| fig4untreated | residual | 0.15 | zero | arithmetic | 0.0085 | 100.0% | 96.2% | 0.1256 | 0.0% | all pass |
| fig4untreated | residual | 0.15 | zero | geometric | 0.0091 | 100.0% | 100.0% | 0.2493 | 0.0% | nominal_coverage_in_bootstrap |
| npc | gaussian | 0.10 | group_common | arithmetic | -0.0007 | 100.0% | 96.2% | 0.0726 | 30.0% | all pass |
| npc | gaussian | 0.10 | group_common | geometric | -0.0005 | 100.0% | 100.0% | 0.1484 | 0.0% | nominal_coverage_in_bootstrap |
| npc | gaussian | 0.15 | group_common | arithmetic | 0.0012 | 100.0% | 92.5% | 0.0888 | 73.8% | all pass |
| npc | gaussian | 0.15 | group_common | geometric | 0.0015 | 100.0% | 100.0% | 0.1721 | 13.8% | nominal_coverage_in_bootstrap |
| npc | gaussian | 0.10 | heterogeneous | arithmetic | -0.0017 | 100.0% | 96.2% | 0.0739 | 70.0% | all pass |
| npc | gaussian | 0.10 | heterogeneous | geometric | -0.0015 | 100.0% | 100.0% | 0.1537 | 3.8% | nominal_coverage_in_bootstrap |
| npc | gaussian | 0.15 | heterogeneous | arithmetic | -0.0038 | 100.0% | 93.8% | 0.0878 | 92.5% | all pass |
| npc | gaussian | 0.15 | heterogeneous | geometric | -0.0035 | 100.0% | 100.0% | 0.1729 | 52.5% | nominal_coverage_in_bootstrap |
| npc | gaussian | 0.10 | homogeneous | arithmetic | 0.0036 | 100.0% | 96.2% | 0.0742 | 25.0% | all pass |
| npc | gaussian | 0.10 | homogeneous | geometric | 0.0037 | 100.0% | 100.0% | 0.1563 | 0.0% | nominal_coverage_in_bootstrap |
| npc | gaussian | 0.15 | homogeneous | arithmetic | -0.0004 | 100.0% | 93.8% | 0.0855 | 78.8% | all pass |
| npc | gaussian | 0.15 | homogeneous | geometric | -0.0002 | 100.0% | 100.0% | 0.1717 | 11.2% | nominal_coverage_in_bootstrap |
| npc | gaussian | 0.10 | zero | arithmetic | 0.0014 | 100.0% | 92.5% | 0.0732 | 0.0% | all pass |
| npc | gaussian | 0.10 | zero | geometric | 0.0014 | 100.0% | 100.0% | 0.1555 | 0.0% | nominal_coverage_in_bootstrap |
| npc | gaussian | 0.15 | zero | arithmetic | 0.0004 | 100.0% | 96.2% | 0.0838 | 0.0% | all pass |
| npc | gaussian | 0.15 | zero | geometric | 0.0006 | 100.0% | 100.0% | 0.1717 | 0.0% | nominal_coverage_in_bootstrap |
| npc | residual | 0.10 | group_common | arithmetic | -0.0020 | 100.0% | 96.2% | 0.0736 | 35.0% | all pass |
| npc | residual | 0.10 | group_common | geometric | -0.0018 | 100.0% | 100.0% | 0.1527 | 1.2% | nominal_coverage_in_bootstrap |
| npc | residual | 0.15 | group_common | arithmetic | -0.0012 | 100.0% | 96.2% | 0.0869 | 78.8% | all pass |
| npc | residual | 0.15 | group_common | geometric | -0.0009 | 100.0% | 100.0% | 0.1691 | 16.2% | nominal_coverage_in_bootstrap |
| npc | residual | 0.10 | heterogeneous | arithmetic | -0.0011 | 100.0% | 95.0% | 0.0741 | 70.0% | all pass |
| npc | residual | 0.10 | heterogeneous | geometric | -0.0010 | 100.0% | 100.0% | 0.1520 | 5.0% | nominal_coverage_in_bootstrap |
| npc | residual | 0.15 | heterogeneous | arithmetic | 0.0000 | 100.0% | 96.2% | 0.0874 | 96.2% | all pass |
| npc | residual | 0.15 | heterogeneous | geometric | 0.0002 | 100.0% | 100.0% | 0.1702 | 52.5% | nominal_coverage_in_bootstrap |
| npc | residual | 0.10 | homogeneous | arithmetic | -0.0001 | 100.0% | 90.0% | 0.0737 | 31.2% | all pass |
| npc | residual | 0.10 | homogeneous | geometric | 0.0001 | 100.0% | 100.0% | 0.1514 | 0.0% | nominal_coverage_in_bootstrap |
| npc | residual | 0.15 | homogeneous | arithmetic | 0.0009 | 100.0% | 93.8% | 0.0860 | 85.0% | all pass |
| npc | residual | 0.15 | homogeneous | geometric | 0.0012 | 100.0% | 100.0% | 0.1675 | 15.0% | nominal_coverage_in_bootstrap |
| npc | residual | 0.10 | zero | arithmetic | -0.0012 | 100.0% | 95.0% | 0.0732 | 0.0% | all pass |
| npc | residual | 0.10 | zero | geometric | -0.0010 | 100.0% | 100.0% | 0.1570 | 0.0% | nominal_coverage_in_bootstrap |
| npc | residual | 0.15 | zero | arithmetic | 0.0016 | 100.0% | 93.8% | 0.0848 | 0.0% | all pass |
| npc | residual | 0.15 | zero | geometric | 0.0018 | 100.0% | 100.0% | 0.1752 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep1 | gaussian | 0.10 | group_common | arithmetic | 0.0007 | 100.0% | 96.2% | 0.0787 | 26.2% | all pass |
| rtx_rep1 | gaussian | 0.10 | group_common | geometric | 0.0009 | 100.0% | 100.0% | 0.1632 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep1 | gaussian | 0.15 | group_common | arithmetic | -0.0113 | 100.0% | 98.8% | 0.1272 | 56.2% | ratio_bias,nominal_coverage_in_bootstrap |
| rtx_rep1 | gaussian | 0.15 | group_common | geometric | -0.0110 | 100.0% | 100.0% | 0.2362 | 5.0% | ratio_bias,nominal_coverage_in_bootstrap |
| rtx_rep1 | gaussian | 0.10 | heterogeneous | arithmetic | -0.0027 | 100.0% | 92.5% | 0.0970 | 87.5% | all pass |
| rtx_rep1 | gaussian | 0.10 | heterogeneous | geometric | -0.0024 | 100.0% | 100.0% | 0.1808 | 47.5% | nominal_coverage_in_bootstrap |
| rtx_rep1 | gaussian | 0.15 | heterogeneous | arithmetic | 0.0030 | 100.0% | 96.2% | 0.1304 | 90.0% | all pass |
| rtx_rep1 | gaussian | 0.15 | heterogeneous | geometric | 0.0038 | 100.0% | 100.0% | 0.2385 | 45.0% | nominal_coverage_in_bootstrap |
| rtx_rep1 | gaussian | 0.10 | homogeneous | arithmetic | -0.0044 | 100.0% | 93.8% | 0.0773 | 40.0% | all pass |
| rtx_rep1 | gaussian | 0.10 | homogeneous | geometric | -0.0041 | 100.0% | 100.0% | 0.1644 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep1 | gaussian | 0.15 | homogeneous | arithmetic | -0.0015 | 100.0% | 95.0% | 0.1182 | 50.0% | all pass |
| rtx_rep1 | gaussian | 0.15 | homogeneous | geometric | -0.0010 | 100.0% | 100.0% | 0.2332 | 6.2% | nominal_coverage_in_bootstrap |
| rtx_rep1 | gaussian | 0.10 | zero | arithmetic | 0.0002 | 100.0% | 93.8% | 0.0807 | 0.0% | all pass |
| rtx_rep1 | gaussian | 0.10 | zero | geometric | 0.0008 | 100.0% | 100.0% | 0.1717 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep1 | gaussian | 0.15 | zero | arithmetic | 0.0010 | 100.0% | 92.5% | 0.1303 | 0.0% | all pass |
| rtx_rep1 | gaussian | 0.15 | zero | geometric | 0.0015 | 100.0% | 100.0% | 0.2544 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep1 | residual | 0.10 | group_common | arithmetic | 0.0004 | 100.0% | 88.8% | 0.0784 | 23.8% | usable_coverage |
| rtx_rep1 | residual | 0.10 | group_common | geometric | 0.0008 | 100.0% | 100.0% | 0.1631 | 1.2% | nominal_coverage_in_bootstrap |
| rtx_rep1 | residual | 0.15 | group_common | arithmetic | 0.0056 | 100.0% | 92.5% | 0.1085 | 46.2% | all pass |
| rtx_rep1 | residual | 0.15 | group_common | geometric | 0.0062 | 100.0% | 100.0% | 0.2111 | 7.5% | nominal_coverage_in_bootstrap |
| rtx_rep1 | residual | 0.10 | heterogeneous | arithmetic | 0.0059 | 100.0% | 91.2% | 0.1052 | 91.2% | all pass |
| rtx_rep1 | residual | 0.10 | heterogeneous | geometric | 0.0065 | 100.0% | 100.0% | 0.1887 | 58.8% | nominal_coverage_in_bootstrap |
| rtx_rep1 | residual | 0.15 | heterogeneous | arithmetic | -0.0003 | 100.0% | 86.2% | 0.1544 | 83.8% | usable_coverage |
| rtx_rep1 | residual | 0.15 | heterogeneous | geometric | 0.0005 | 100.0% | 98.8% | 0.2543 | 82.5% | nominal_coverage_in_bootstrap |
| rtx_rep1 | residual | 0.10 | homogeneous | arithmetic | -0.0010 | 100.0% | 96.2% | 0.0794 | 31.2% | all pass |
| rtx_rep1 | residual | 0.10 | homogeneous | geometric | -0.0009 | 100.0% | 100.0% | 0.1630 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep1 | residual | 0.15 | homogeneous | arithmetic | 0.0016 | 100.0% | 97.5% | 0.1165 | 46.2% | all pass |
| rtx_rep1 | residual | 0.15 | homogeneous | geometric | 0.0024 | 100.0% | 100.0% | 0.2178 | 3.8% | nominal_coverage_in_bootstrap |
| rtx_rep1 | residual | 0.10 | zero | arithmetic | -0.0021 | 100.0% | 92.5% | 0.0724 | 0.0% | all pass |
| rtx_rep1 | residual | 0.10 | zero | geometric | -0.0016 | 100.0% | 100.0% | 0.1585 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep1 | residual | 0.15 | zero | arithmetic | -0.0023 | 100.0% | 98.8% | 0.1028 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep1 | residual | 0.15 | zero | geometric | -0.0014 | 100.0% | 100.0% | 0.2116 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep2 | gaussian | 0.10 | group_common | arithmetic | 0.0005 | 100.0% | 93.8% | 0.0495 | 56.2% | all pass |
| rtx_rep2 | gaussian | 0.10 | group_common | geometric | 0.0008 | 100.0% | 100.0% | 0.1114 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep2 | gaussian | 0.15 | group_common | arithmetic | 0.0025 | 100.0% | 95.0% | 0.0740 | 78.8% | all pass |
| rtx_rep2 | gaussian | 0.15 | group_common | geometric | 0.0030 | 100.0% | 100.0% | 0.1573 | 26.2% | nominal_coverage_in_bootstrap |
| rtx_rep2 | gaussian | 0.10 | heterogeneous | arithmetic | -0.0079 | 100.0% | 92.5% | 0.0784 | 92.5% | all pass |
| rtx_rep2 | gaussian | 0.10 | heterogeneous | geometric | -0.0078 | 100.0% | 100.0% | 0.1336 | 97.5% | nominal_coverage_in_bootstrap |
| rtx_rep2 | gaussian | 0.15 | heterogeneous | arithmetic | -0.0019 | 100.0% | 93.8% | 0.0966 | 92.5% | all pass |
| rtx_rep2 | gaussian | 0.15 | heterogeneous | geometric | -0.0017 | 100.0% | 100.0% | 0.1574 | 98.8% | nominal_coverage_in_bootstrap |
| rtx_rep2 | gaussian | 0.10 | homogeneous | arithmetic | 0.0002 | 100.0% | 95.0% | 0.0516 | 60.0% | all pass |
| rtx_rep2 | gaussian | 0.10 | homogeneous | geometric | 0.0005 | 100.0% | 100.0% | 0.1193 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep2 | gaussian | 0.15 | homogeneous | arithmetic | 0.0008 | 100.0% | 100.0% | 0.0784 | 88.8% | nominal_coverage_in_bootstrap |
| rtx_rep2 | gaussian | 0.15 | homogeneous | geometric | 0.0012 | 100.0% | 100.0% | 0.1590 | 18.8% | nominal_coverage_in_bootstrap |
| rtx_rep2 | gaussian | 0.10 | zero | arithmetic | -0.0023 | 100.0% | 96.2% | 0.0482 | 0.0% | all pass |
| rtx_rep2 | gaussian | 0.10 | zero | geometric | -0.0021 | 100.0% | 100.0% | 0.1139 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep2 | gaussian | 0.15 | zero | arithmetic | -0.0010 | 100.0% | 96.2% | 0.0703 | 0.0% | all pass |
| rtx_rep2 | gaussian | 0.15 | zero | geometric | -0.0008 | 100.0% | 100.0% | 0.1483 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep2 | residual | 0.10 | group_common | arithmetic | -0.0034 | 100.0% | 96.2% | 0.0501 | 67.5% | all pass |
| rtx_rep2 | residual | 0.10 | group_common | geometric | -0.0032 | 100.0% | 100.0% | 0.1140 | 3.8% | nominal_coverage_in_bootstrap |
| rtx_rep2 | residual | 0.15 | group_common | arithmetic | -0.0034 | 100.0% | 93.8% | 0.0711 | 86.2% | all pass |
| rtx_rep2 | residual | 0.15 | group_common | geometric | -0.0032 | 100.0% | 100.0% | 0.1431 | 31.2% | nominal_coverage_in_bootstrap |
| rtx_rep2 | residual | 0.10 | heterogeneous | arithmetic | -0.0042 | 100.0% | 93.8% | 0.1087 | 93.8% | all pass |
| rtx_rep2 | residual | 0.10 | heterogeneous | geometric | -0.0039 | 100.0% | 98.8% | 0.1587 | 96.2% | nominal_coverage_in_bootstrap |
| rtx_rep2 | residual | 0.15 | heterogeneous | arithmetic | 0.0029 | 100.0% | 91.2% | 0.1230 | 91.2% | all pass |
| rtx_rep2 | residual | 0.15 | heterogeneous | geometric | 0.0032 | 100.0% | 98.8% | 0.1824 | 98.8% | nominal_coverage_in_bootstrap |
| rtx_rep2 | residual | 0.10 | homogeneous | arithmetic | 0.0023 | 100.0% | 95.0% | 0.0473 | 63.7% | all pass |
| rtx_rep2 | residual | 0.10 | homogeneous | geometric | 0.0025 | 100.0% | 100.0% | 0.1086 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep2 | residual | 0.15 | homogeneous | arithmetic | -0.0012 | 100.0% | 95.0% | 0.0701 | 86.2% | all pass |
| rtx_rep2 | residual | 0.15 | homogeneous | geometric | -0.0010 | 100.0% | 100.0% | 0.1447 | 37.5% | nominal_coverage_in_bootstrap |
| rtx_rep2 | residual | 0.10 | zero | arithmetic | -0.0018 | 100.0% | 96.2% | 0.0449 | 0.0% | all pass |
| rtx_rep2 | residual | 0.10 | zero | geometric | -0.0017 | 100.0% | 100.0% | 0.1108 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep2 | residual | 0.15 | zero | arithmetic | -0.0004 | 100.0% | 93.8% | 0.0642 | 0.0% | all pass |
| rtx_rep2 | residual | 0.15 | zero | geometric | -0.0002 | 100.0% | 100.0% | 0.1361 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | gaussian | 0.10 | group_common | arithmetic | 0.0019 | 100.0% | 95.0% | 0.0679 | 32.5% | all pass |
| rtx_rep3 | gaussian | 0.10 | group_common | geometric | 0.0020 | 100.0% | 100.0% | 0.1463 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | gaussian | 0.15 | group_common | arithmetic | 0.0024 | 100.0% | 95.0% | 0.1116 | 50.0% | all pass |
| rtx_rep3 | gaussian | 0.15 | group_common | geometric | 0.0028 | 100.0% | 100.0% | 0.2115 | 5.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | gaussian | 0.10 | heterogeneous | arithmetic | 0.0011 | 100.0% | 97.5% | 0.0728 | 92.5% | all pass |
| rtx_rep3 | gaussian | 0.10 | heterogeneous | geometric | 0.0013 | 100.0% | 100.0% | 0.1438 | 32.5% | nominal_coverage_in_bootstrap |
| rtx_rep3 | gaussian | 0.15 | heterogeneous | arithmetic | -0.0045 | 100.0% | 91.2% | 0.1063 | 88.8% | all pass |
| rtx_rep3 | gaussian | 0.15 | heterogeneous | geometric | -0.0042 | 100.0% | 100.0% | 0.1990 | 60.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | gaussian | 0.10 | homogeneous | arithmetic | 0.0034 | 100.0% | 95.0% | 0.0699 | 30.0% | all pass |
| rtx_rep3 | gaussian | 0.10 | homogeneous | geometric | 0.0035 | 100.0% | 100.0% | 0.1535 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | gaussian | 0.15 | homogeneous | arithmetic | 0.0029 | 100.0% | 92.5% | 0.1075 | 52.5% | all pass |
| rtx_rep3 | gaussian | 0.15 | homogeneous | geometric | 0.0033 | 100.0% | 100.0% | 0.2077 | 10.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | gaussian | 0.10 | zero | arithmetic | -0.0004 | 100.0% | 88.8% | 0.0676 | 0.0% | usable_coverage |
| rtx_rep3 | gaussian | 0.10 | zero | geometric | -0.0001 | 100.0% | 100.0% | 0.1506 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | gaussian | 0.15 | zero | arithmetic | 0.0053 | 100.0% | 95.0% | 0.1024 | 0.0% | all pass |
| rtx_rep3 | gaussian | 0.15 | zero | geometric | 0.0058 | 100.0% | 100.0% | 0.2147 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | residual | 0.10 | group_common | arithmetic | -0.0026 | 100.0% | 95.0% | 0.0680 | 38.8% | all pass |
| rtx_rep3 | residual | 0.10 | group_common | geometric | -0.0022 | 100.0% | 100.0% | 0.1446 | 2.5% | nominal_coverage_in_bootstrap |
| rtx_rep3 | residual | 0.15 | group_common | arithmetic | 0.0011 | 100.0% | 93.8% | 0.0995 | 66.2% | all pass |
| rtx_rep3 | residual | 0.15 | group_common | geometric | 0.0016 | 100.0% | 100.0% | 0.1871 | 7.5% | nominal_coverage_in_bootstrap |
| rtx_rep3 | residual | 0.10 | heterogeneous | arithmetic | 0.0027 | 100.0% | 95.0% | 0.0790 | 93.8% | all pass |
| rtx_rep3 | residual | 0.10 | heterogeneous | geometric | 0.0029 | 100.0% | 100.0% | 0.1545 | 46.2% | nominal_coverage_in_bootstrap |
| rtx_rep3 | residual | 0.15 | heterogeneous | arithmetic | 0.0087 | 100.0% | 90.0% | 0.1195 | 88.8% | all pass |
| rtx_rep3 | residual | 0.15 | heterogeneous | geometric | 0.0092 | 100.0% | 100.0% | 0.2049 | 82.5% | nominal_coverage_in_bootstrap |
| rtx_rep3 | residual | 0.10 | homogeneous | arithmetic | 0.0010 | 100.0% | 95.0% | 0.0652 | 32.5% | all pass |
| rtx_rep3 | residual | 0.10 | homogeneous | geometric | 0.0013 | 100.0% | 100.0% | 0.1499 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | residual | 0.15 | homogeneous | arithmetic | 0.0025 | 100.0% | 92.5% | 0.0956 | 65.0% | all pass |
| rtx_rep3 | residual | 0.15 | homogeneous | geometric | 0.0030 | 100.0% | 100.0% | 0.1890 | 15.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | residual | 0.10 | zero | arithmetic | 0.0005 | 100.0% | 93.8% | 0.0650 | 0.0% | all pass |
| rtx_rep3 | residual | 0.10 | zero | geometric | 0.0007 | 100.0% | 100.0% | 0.1456 | 0.0% | nominal_coverage_in_bootstrap |
| rtx_rep3 | residual | 0.15 | zero | arithmetic | 0.0056 | 100.0% | 96.2% | 0.0897 | 0.0% | all pass |
| rtx_rep3 | residual | 0.15 | zero | geometric | 0.0060 | 100.0% | 100.0% | 0.1983 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | gaussian | 0.10 | group_common | arithmetic | 0.0015 | 100.0% | 95.0% | 0.0765 | 31.2% | all pass |
| sdc_corrected | gaussian | 0.10 | group_common | geometric | 0.0010 | 100.0% | 100.0% | 0.1734 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | gaussian | 0.15 | group_common | arithmetic | -0.0004 | 100.0% | 97.5% | 0.2013 | 16.2% | all pass |
| sdc_corrected | gaussian | 0.15 | group_common | geometric | 0.0026 | 100.0% | 100.0% | 0.3849 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | gaussian | 0.10 | heterogeneous | arithmetic | -0.0022 | 100.0% | 97.5% | 0.0805 | 97.5% | all pass |
| sdc_corrected | gaussian | 0.10 | heterogeneous | geometric | -0.0024 | 100.0% | 100.0% | 0.1729 | 61.3% | nominal_coverage_in_bootstrap |
| sdc_corrected | gaussian | 0.15 | heterogeneous | arithmetic | -0.0086 | 100.0% | 97.5% | 0.2012 | 75.0% | all pass |
| sdc_corrected | gaussian | 0.15 | heterogeneous | geometric | -0.0076 | 100.0% | 100.0% | 0.3730 | 8.8% | nominal_coverage_in_bootstrap |
| sdc_corrected | gaussian | 0.10 | homogeneous | arithmetic | -0.0049 | 100.0% | 96.2% | 0.0753 | 35.0% | all pass |
| sdc_corrected | gaussian | 0.10 | homogeneous | geometric | -0.0038 | 100.0% | 100.0% | 0.1706 | 1.2% | nominal_coverage_in_bootstrap |
| sdc_corrected | gaussian | 0.15 | homogeneous | arithmetic | -0.0061 | 100.0% | 87.5% | 0.1999 | 21.2% | usable_coverage |
| sdc_corrected | gaussian | 0.15 | homogeneous | geometric | -0.0049 | 100.0% | 100.0% | 0.3705 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | gaussian | 0.10 | zero | arithmetic | 0.0044 | 100.0% | 96.2% | 0.0757 | 0.0% | all pass |
| sdc_corrected | gaussian | 0.10 | zero | geometric | 0.0051 | 100.0% | 100.0% | 0.1769 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | gaussian | 0.15 | zero | arithmetic | -0.0095 | 100.0% | 98.8% | 0.2093 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | gaussian | 0.15 | zero | geometric | -0.0064 | 100.0% | 100.0% | 0.4008 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | residual | 0.10 | group_common | arithmetic | 0.0000 | 100.0% | 93.8% | 0.0264 | 93.8% | all pass |
| sdc_corrected | residual | 0.10 | group_common | geometric | 0.0003 | 100.0% | 100.0% | 0.0654 | 31.2% | nominal_coverage_in_bootstrap |
| sdc_corrected | residual | 0.15 | group_common | arithmetic | -0.0023 | 100.0% | 100.0% | 0.0417 | 100.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | residual | 0.15 | group_common | geometric | -0.0021 | 100.0% | 100.0% | 0.0877 | 98.8% | nominal_coverage_in_bootstrap |
| sdc_corrected | residual | 0.10 | heterogeneous | arithmetic | 0.0016 | 100.0% | 95.0% | 0.0506 | 95.0% | all pass |
| sdc_corrected | residual | 0.10 | heterogeneous | geometric | 0.0015 | 100.0% | 100.0% | 0.0807 | 100.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | residual | 0.15 | heterogeneous | arithmetic | -0.0007 | 100.0% | 91.2% | 0.0831 | 91.2% | all pass |
| sdc_corrected | residual | 0.15 | heterogeneous | geometric | -0.0004 | 100.0% | 97.5% | 0.1135 | 97.5% | all pass |
| sdc_corrected | residual | 0.10 | homogeneous | arithmetic | 0.0007 | 100.0% | 93.8% | 0.0272 | 93.8% | all pass |
| sdc_corrected | residual | 0.10 | homogeneous | geometric | 0.0011 | 100.0% | 100.0% | 0.0671 | 22.5% | nominal_coverage_in_bootstrap |
| sdc_corrected | residual | 0.15 | homogeneous | arithmetic | -0.0002 | 100.0% | 92.5% | 0.0419 | 92.5% | all pass |
| sdc_corrected | residual | 0.15 | homogeneous | geometric | 0.0005 | 100.0% | 100.0% | 0.0915 | 95.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | residual | 0.10 | zero | arithmetic | 0.0005 | 100.0% | 97.5% | 0.0250 | 0.0% | all pass |
| sdc_corrected | residual | 0.10 | zero | geometric | 0.0003 | 100.0% | 100.0% | 0.0650 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_corrected | residual | 0.15 | zero | arithmetic | 0.0002 | 100.0% | 95.0% | 0.0370 | 0.0% | all pass |
| sdc_corrected | residual | 0.15 | zero | geometric | -0.0000 | 100.0% | 100.0% | 0.0877 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_raw | gaussian | 0.10 | group_common | arithmetic | 0.0027 | 100.0% | 96.2% | 0.0802 | 18.8% | all pass |
| sdc_raw | gaussian | 0.10 | group_common | geometric | 0.0032 | 100.0% | 100.0% | 0.1804 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_raw | gaussian | 0.15 | group_common | arithmetic | 0.0074 | 100.0% | 96.2% | 0.2256 | 6.2% | all pass |
| sdc_raw | gaussian | 0.15 | group_common | geometric | 0.0092 | 100.0% | 100.0% | 0.4131 | 1.2% | nominal_coverage_in_bootstrap |
| sdc_raw | gaussian | 0.10 | heterogeneous | arithmetic | 0.0054 | 100.0% | 95.0% | 0.0827 | 95.0% | all pass |
| sdc_raw | gaussian | 0.10 | heterogeneous | geometric | 0.0063 | 100.0% | 100.0% | 0.1750 | 33.8% | nominal_coverage_in_bootstrap |
| sdc_raw | gaussian | 0.15 | heterogeneous | arithmetic | 0.0041 | 100.0% | 96.2% | 0.2246 | 52.5% | all pass |
| sdc_raw | gaussian | 0.15 | heterogeneous | geometric | 0.0068 | 100.0% | 100.0% | 0.4019 | 3.8% | nominal_coverage_in_bootstrap |
| sdc_raw | gaussian | 0.10 | homogeneous | arithmetic | -0.0001 | 100.0% | 95.0% | 0.0785 | 27.5% | all pass |
| sdc_raw | gaussian | 0.10 | homogeneous | geometric | 0.0002 | 100.0% | 100.0% | 0.1765 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_raw | gaussian | 0.15 | homogeneous | arithmetic | 0.0080 | 100.0% | 95.0% | 0.2225 | 17.5% | all pass |
| sdc_raw | gaussian | 0.15 | homogeneous | geometric | 0.0118 | 100.0% | 100.0% | 0.4180 | 0.0% | ratio_bias,nominal_coverage_in_bootstrap |
| sdc_raw | gaussian | 0.10 | zero | arithmetic | -0.0013 | 100.0% | 93.8% | 0.0777 | 0.0% | all pass |
| sdc_raw | gaussian | 0.10 | zero | geometric | -0.0013 | 100.0% | 100.0% | 0.1788 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_raw | gaussian | 0.15 | zero | arithmetic | -0.0003 | 100.0% | 96.2% | 0.2304 | 0.0% | all pass |
| sdc_raw | gaussian | 0.15 | zero | geometric | 0.0032 | 100.0% | 100.0% | 0.4337 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_raw | residual | 0.10 | group_common | arithmetic | -0.0004 | 100.0% | 96.2% | 0.0271 | 96.2% | all pass |
| sdc_raw | residual | 0.10 | group_common | geometric | -0.0008 | 100.0% | 100.0% | 0.0653 | 31.2% | nominal_coverage_in_bootstrap |
| sdc_raw | residual | 0.15 | group_common | arithmetic | -0.0013 | 100.0% | 93.8% | 0.0424 | 93.8% | all pass |
| sdc_raw | residual | 0.15 | group_common | geometric | -0.0022 | 100.0% | 100.0% | 0.0887 | 93.8% | nominal_coverage_in_bootstrap |
| sdc_raw | residual | 0.10 | heterogeneous | arithmetic | 0.0009 | 100.0% | 93.8% | 0.0496 | 93.8% | all pass |
| sdc_raw | residual | 0.10 | heterogeneous | geometric | 0.0011 | 100.0% | 100.0% | 0.0797 | 100.0% | nominal_coverage_in_bootstrap |
| sdc_raw | residual | 0.15 | heterogeneous | arithmetic | 0.0048 | 100.0% | 93.8% | 0.0802 | 93.8% | all pass |
| sdc_raw | residual | 0.15 | heterogeneous | geometric | 0.0050 | 100.0% | 97.5% | 0.1107 | 97.5% | all pass |
| sdc_raw | residual | 0.10 | homogeneous | arithmetic | -0.0000 | 100.0% | 95.0% | 0.0272 | 95.0% | all pass |
| sdc_raw | residual | 0.10 | homogeneous | geometric | -0.0001 | 100.0% | 100.0% | 0.0660 | 23.8% | nominal_coverage_in_bootstrap |
| sdc_raw | residual | 0.15 | homogeneous | arithmetic | 0.0031 | 100.0% | 88.8% | 0.0423 | 88.8% | usable_coverage |
| sdc_raw | residual | 0.15 | homogeneous | geometric | 0.0031 | 100.0% | 100.0% | 0.0905 | 90.0% | nominal_coverage_in_bootstrap |
| sdc_raw | residual | 0.10 | zero | arithmetic | -0.0007 | 100.0% | 96.2% | 0.0252 | 0.0% | all pass |
| sdc_raw | residual | 0.10 | zero | geometric | -0.0004 | 100.0% | 100.0% | 0.0651 | 0.0% | nominal_coverage_in_bootstrap |
| sdc_raw | residual | 0.15 | zero | arithmetic | -0.0008 | 100.0% | 96.2% | 0.0378 | 0.0% | all pass |
| sdc_raw | residual | 0.15 | zero | geometric | -0.0003 | 100.0% | 100.0% | 0.0881 | 0.0% | nominal_coverage_in_bootstrap |
