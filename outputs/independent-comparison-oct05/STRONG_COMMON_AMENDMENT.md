# Information-matched common-round comparator before validation

Exposed development now includes the unperturbed crop and common-round doses
1 and 2 (both seeds), plus dose 4 seed 501. The original phase-only common-round
policy can lose about .07–.09 identity F1 even when the added shift is truly
common across all four colors. Comparing per-frame iteration only against that
policy could reward its optimizer, rather than its extra calibration degrees.
No color-offset condition and no validation crop outcome has been read at this
amendment. Retain the initial policies and add this stronger comparator.

Add `common_iterative`: sum four fixed-preprocessed channels in each round,
initialize from the original common-round phase shifts, and perform ten fixed
population-template updates across the seven round sums. Use the same bounded
amplitude correlation and original-frame resampling as the per-frame template
policy. Center round shifts and repeat each correction over its four channels.
Use the same detection thresholds and endpoints; all its failures and runtime
count. This matches the optimizer and inputs while restricting corrections to
common round translations. The per-frame policy stays unchanged.

The original observable decision extension remains reported against phase-only
alignment as specified. A separate **strong-comparator** rule uses RMS difference
between per-frame and common-iterative corrections, the same fixed .05 F1 target,
threshold grid and development qualification gate. Both rules must freeze before
the validation crop is opened. The main practical interpretation must use the
strong-comparator result; improvement over phase-only cannot substitute for it.

This is six policies, twelve decoder arms per condition, 312 final arms across
the two crops. Initial partial development template failures are additional
retained evidence, not among these arms. There is still only one documented
BARseq benchmark source preparation and no population replication claim.
