# Matched-distribution public recovery, before perturbed decode outcomes

Main question: does the placement of geometric error across rounds versus
colors change molecular-output recovery, even when the complete marginal
offset distribution, two-dimensional covariance and RMS are identical?
Strong alternatives may match
the per-frame method. The comparative conclusion does not require a unique
method win. Also test the existing >=2-pixel image-disagreement recommendation,
without fitting its cut on this preparation.

This comparison runs only if both fixed native fields pass both decoder gates
and the preceding integer-grid control verifies exact oracle calls. The prior
partial adapter and its failed gates remain separate. Source zero is exposed
adapter development; source one's new fixed crop is a spatial holdout. Both
are one BARseq3 preparation, distinct from the older CSHL annotated benchmark;
they are not two biological preparations or independent prospective validation.

Fixed seed 505, ten rounds/four colors, vector v=(4,2). Common-only displacement
assigns +v to five rounds and -v to five rounds, with all four colors sharing
each round's offset. Color-only displacement assigns +v to two colors and -v
to two colors across all rounds. Shuffle round order and color assignment
using the seed before any outcome. Both complete 40-frame arrays contain
exactly twenty +v and twenty -v: zero mean, covariance [[16,8],[8,4]], RMS
sqrt(20), identical tails, maxima and every marginal moment. The common versus
color grouping alone changes. Covariance is intentionally rank one; this
minimal test does not establish generality over arbitrary displacement fields.
Added offsets and every applied fitted correction are integer pixels.

Pre-execution source review replaces the unrun covariance-only design in
commit a8e3c42. That earlier pair matched covariance/RMS but had different
marginal offset distributions. The stronger design removes that avoidable
confound before any recovery fit or decode outcome. No result is discarded.

Policies: none; ten-step strong common-round template; ten-step per-frame
template with the same optimizer; actual author barcode-model optimizer,
50 steps from two image-only initializations selected by its final objective;
inverse known added offsets. Fit on the same 128-pixel prepared native input.
The native conditioning is performed before adding offsets. No second global
normalization or unmixing is applied. Decoder kernels remain the stock
correlation (.775 cut, .05 noise floor) and the actual 100-step nonnegative
regression (.20 of maximum, zero PSF/L1, estimated gains).

Decode a common [16,112) window, restore its full crop coordinates, and score
[24,104). The zero/no-correction window supplies a decoder-specific reference
and must pass the unchanged >=50 nonblank / <=.10 proxy / >=99% support gates.
It is an algorithm-output reference, not independent RNA identity truth. The
oracle must reproduce native calls exactly; the shared decoder domain controls
gain and maximum-threshold changes from discarded outer pixels.
Oracle and zero/no-correction arms reuse the reference kernel only after exact
equality of the full decoder input window is asserted. They are retained arms,
not independent solver executions or additional biological replication. The
measured reference kernel cost remains identified separately from execution cost.
Spatial matching is label-independent maximum-cardinality/minimum-distance at 5 pixels,
with 1.5/3-pixel and label-constrained sensitivities. Retain gene-class identity
F1, composition TV, yield, blank proxies, map residuals and all support failures.

Two fields x three conditions (zero, common, round-plus-color) x five policies
x two kernels = 60 planned arms. Count every missing, failed and support-rejected
arm. Report each field/family separately; full-cohort conclusions require all
60 results. Primary benefit is >=.05 identity-F1 versus strong common-round;
harm is <=-.05. The existing 2-pixel recommendation uses only fitted-map
disagreement, never imposed truth. Runtime includes both fits needed for that
feature, both model initializations, preprocessing/warp/kernel/peaks, and is
descriptive on local hardware. No speed or true-FDR superiority is inferred.

Do not tune a threshold, choose another crop, drop a decoder, soften a gate or
change the perturbation after outcomes. A failed or negative comparison remains
a result. Even a success supports retrospective controlled recovery of native
outputs, with annotation-grounded identity evaluated separately on CSHL data.
It cannot validate the physical randomized temporal intervention, resolve the
failed S3/ISS confirmations, or guarantee Nature Methods acceptance.
