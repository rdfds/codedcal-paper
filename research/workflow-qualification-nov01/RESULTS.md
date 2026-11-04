# Completed development comparison: no demonstrated readout-aware advantage

The frozen matrix is complete. It does **not** support the proposed improvement
in workflow qualification. On the ten workflows with supported held-out outcomes,
the mixed conventional suite detects all nine failures and accepts the sole
adequate workflow. The readout-aware suite detects five failures, accepts four
failed workflows and rejects the adequate workflow. No fresh comparative outcomes
are authorized by this result. No estimator or manuscript claim is upgraded.

| Four-pair probing suite | Failed workflows rejected / 9 | Adequate workflows accepted / 1 |
| --- | ---: | ---: |
| Conventional round/color | 5 | 1 |
| Mixed conventional | 9 | 1 |
| Balanced random | 9 | 0 |
| Readout-aware coverage heuristic | 5 | 0 |

These are descriptive cases, not independent replicates or estimates of practical
error prevalence. There is only one adequate example, below the frozen minimum
of five adequate examples. Leave-preparation-out scalar QC and combined QC rules
therefore abstain on all ten supported cases; that is not evidence that probing
outperforms a properly trained QC rule. The digital held-out bank has five known
challenges from the planned families. Passing it does not establish safety under
unknown natural errors or all possible transformations.

## Complete accounting

All 32 original workflows are retained: eight image-regression compute exclusions,
eight source-gate exclusions, and sixteen completed calibration/decoder workflows.
Six of those sixteen fail the complete native workflow gate and have unsupported
held-out truth. The ten supported outcomes arise from two preparations: two fields
of the BARseq3 acquisition, and the conditional STARfinder tissue target with two
decoders. Fields, decoders and challenges do not add independent preparations.

All 332 challenge-pair records and sixteen native records are retained unchanged.
The full union dispatched 680 logical workflow queries: 380 executed workflows
and 300 exact-result cache hits. They performed 284 actual calibration executions
and 380 actual decoder executions. The frozen maximum reservation is 752 logical
calls, including strategy-overlapping probes. Do not confuse that reservation
with dispatched or actual executions.

| Calibration policy | Actual workflows | Actual calibrations | Actual decoders | Exact-result hits | Measured query-union seconds |
| --- | ---: | ---: | ---: | ---: | ---: |
| Round | 95 | 71 | 95 | 75 | 102.49 |
| Frame template | 95 | 71 | 95 | 75 | 824.75 |
| Barcode model | 95 | 71 | 95 | 75 | 4207.09 |
| Group search | 95 | 71 | 95 | 75 | 1186.66 |

These are measured unions of all queries, including shared cache lookup work;
they are not separately measured per-strategy times. Maximum invocation budgets
are matched, not observed wall-clock cost. The worker records 6659.74 elapsed
monotonic seconds; host sleep separates that counter from the UTC launch-to-exit
span. Neither counter is CPU time. No paid or cloud work was performed.

## Reference audit and retained failure evidence

Every completed integer oracle catalog and correction map exactly equals its
corresponding complete native workflow. Reference-first offline adjudication
changes zero status/failure classifications. Raw v1 assessments are retained in
the original JSON rows; subsequent software repairs do not refit those outcomes.
Unresolved spatial-first matching is unavailable, rather than a computed score.
The label-aware primary matching and composition endpoint remain separate.

`development-raw-manifest.json` seals the raw rows and log. `frozen-v1/` retains
the exact sources pinned by `development-freeze.json`, including the older API,
so its implementation is not lost when the shipped API is corrected. The frozen
run is not restarted or replaced. `adjudicated-summary.json`,
`workflow-decisions.csv` and `recovery-profiles.csv` provide the complete analysis.

The separate public culture native qualification is bounded by
`FULL_WORKFLOW_NATIVE_QUALIFICATION.md` and `NATIVE_SAMPLE_SUPPORT.md`. It is not
a comparative confirmation and cannot reverse the failed development utility
gate. Its earlier 256-pixel failures remain evidence. No physical-motion,
independent RNA-truth or randomized molecular-readout validation follows here.

## Bounded native qualification and software completion

The one fixed 512-pixel hybrid workflow completed all four native control calls
in 141.17 recorded monotonic seconds: four logical calls, two actual workflows
and two exact-result cache hits. Native/no-effect/integer-inverse catalogs each
contain 857 coding calls and zero named blanks, with output support 1.0. Their
native blank-proxy upper95 is 0.07806848, meeting the frozen 0.10 bound. Integer
pixels/support and native/no-effect/inverse catalogs are exact. The portable
processed input, original acquisition receipts, actual source pins and all
per-call checkpoints are retained. This is an adapted existing-kernel hybrid,
not complete published Serval reproduction or independently established accuracy.

The composed fractional control has 661 coding calls and zero blanks. Its proxy
upper95 is approximately 0.1012: it has enough counts/support for the computational
reference floor, but does **not** independently meet the native 0.10 quality
bound. The raw script's `eligible_for_confirmation` flag describes its native
gate, not successful development utility or fractional molecular quality. The
separate `confirmation-decision.json` explicitly authorizes no comparative calls.
There is no additional crop, method or parameter search.

The shipped new API now checks reference counts before candidate execution/
support failures, marks unresolved spatial-first matching unavailable, and accepts
an explicit complete native policy. The separate workflow CLI writes readable
recovery/composition/support/failure reports and qualifies each scored policy's
native catalog. The original estimator CLI remains byte-identical. All 261
package/adjudication checks pass. These repairs do not change a fitted catalog,
the preserved raw v1 assessments, the comparative classifications or runtimes.

The conservative total charge is 843 of 1600: 84 earlier control calls, two
calibration-only and one decoder-only synthetic controls, the complete 752-call
development reservation, and four full-workflow native controls. There are zero
external comparative calls and zero new network requests in this completion.
