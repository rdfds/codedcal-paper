# Focused reconciliation sensitivity — specified 23 October 2026

Post-hoc sensitivity of exposed Source Data 69, not independent confirmation.
This specification is committed before new sensitivity outcomes. Historical
packages, cohort, raw filters, event linking, temporal panels, frequencies,
code assignments and previous failures are immutable. No D076 data are used.

## Fixed cohort and grid

Include all eight acquisitions / nine views in the original registry: NPC,
EDF10 replicates 1–2, Figure4 untreated, RTX replicates 1–3, and SDC corrected
and raw. SDC views are paired observations of one acquisition; Figure4 RTX is
an alias of RTX replicate 1, not an additional acquisition.

Cross three strict single-linkage radii (0.8, 1.0, 1.2 times the historical
10 nm NPC/CD20 or 8 nm SDC radius) with three eligibility profiles:

| Parameter | Lenient | Original | Stringent |
|---|---:|---:|---:|
| Minimum events | ceiling(0.8 × original) | original | ceiling(1.2 × original) |
| Minimum events in each panel | original − 1, at least 1 | original | original + 1 |
| Maximum centered RMS | 1.2 × original | original | 0.8 × original |
| Minimum source separation | 0.8 × original | original | 1.2 × original |
| Maximum split gain | original + 0.05, at most 1 | original | original − 0.05 |
| Maximum z spread, if applicable | 1.2 × original | original | 0.8 × original |

The nine crossed settings per view are a bounded neighborhood, not confidence
bounds on source identity or a validated set of physical alternatives. No
setting is selected for agreement. Binding-event formation and localization
precision gates remain fixed; their sensitivity is outside this experiment.
Report three quadratic reference weightings (uniform, event-count,
localization-count) in every arm. Preserve the distinct finite-frequency
diagonal target; do not interpret quadratic-versus-Fourier differences as bias.

## Endpoints and uncertainty

Reconstruct source-to-localization membership using the existing event mapping.
Replay the original setting against all retained per-group code and diagonal
arrays, source counts, weights and quadratic estimates before interpreting
other settings. Compare partitions rather than arbitrary cluster label numbers.
Use centered float64 throughout new arms. Separate the historical-to-float64
change from float64 full-population-to-selected changes.

For each setting retain source counts, rejection reasons, event/localization
coverage, overlap with original selected rows, effective weights, selected code,
matched diagonal and full-population moments. The descriptive fit-gap question
is whether restricting code to the same rows reduces its absolute gap from the
diagonal reference. Retain worsening, unavailable and optimizer-boundary cases.
Do not promote a bounded two-frequency fit into a physical uncertainty interval.

At both original frequencies (0.10, 0.15 nm^-1), retain Fieller ratio sets,
including disconnected, unbounded and empty outcomes. Transform the entire
positive-ratio portion to signed equivalent-T sets without clipping to positive
T. Use the original fixed 4 × 3 spatial partition, paired 2,000 block-bootstrap
draws and leave-one-block-out estimates. Preserve invalid draws; report a
paired percentile interval only with at least 95% valid draws and mark it
conditional/pointwise. Spatial blocks are not independent preparations.

Summarize all nine settings per view, the range of gap changes, residual
frequency-specific disagreements, weight sensitivity, and every unsupported
case. A grid majority is descriptive, not a replication p-value. A zero-spanning
interval is not equivalence. Recount the original 7/9 historical comparison
separately from the numerically standardized sensitivity results.

## Existing response-fit uncertainty

The six CD20 preparations have retained paired per-group control arrays for
two noise laws, four imposed T values and two seeds. Reconstruct all fitted
intercepts/slopes and their covariance, preserving within-seed pairing across
T. Propagate observation and response-coefficient uncertainty jointly using a
Bonferroni rectangle: observation confidence 97.5%, intercept and slope each
98.75%, yielding nominal joint coverage at least 95% under their respective
approximations without assuming independence. Transform all rectangle corners
when slope is strictly positive; retain an unsupported/unbounded result if it
crosses zero. Include each noise model separately and their sensitivity envelope.
Retain the original conditional qualification gates and failures.

These intervals include finite observation and Monte Carlo fit uncertainty,
conditional on inferred templates. They do not include unknown physical source
assignments, all template uncertainty or all possible noise laws. NPC and SDC
do not have these same preparation-specific response controls; do not transfer
CD20 fits or fabricate calibrated intervals for them. All nine views receive
the full observed-data grid and uncertainty analysis above. Response fits are
not rerun under altered QC, so calibrated response sensitivities apply only to
the original reference, never to every new grid arm.

## Verification and reporting

One computationally heavy process at a time. No paid jobs or new downloads.
Hash every original dependency and verify it before use; record software and
run status. Preserve per-group aggregates sufficient for independent recounts.
Test strict-cutoff partitions, boundary-crossing event membership, empty source
sets, ratio-set topology and joint rectangle transformations. Independently
recount all 81 arms and baseline replays. Correct implementation failures in
new code without modifying old evidence or this grid to favor agreement.

Update the manuscript and source ledger to the scope actually supported. Keep
original negative results and theoretical limits. Small commits and pushes;
no merge or submission. Computational robustness cannot establish physical
temporal truth, validate unselected sources or guarantee journal acceptance.
