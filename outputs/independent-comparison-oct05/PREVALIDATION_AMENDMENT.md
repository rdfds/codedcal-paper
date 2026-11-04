# Strong-alternative amendment, before perturbed or validation outcomes

The unperturbed development qualification passed for both families (identity
F1 .8643 correlation and .8571 regression). Their selected thresholds remain
.85 absolute and .20 of the density maximum, respectively. No nonzero condition
and no validation image/annotation value has been analyzed at this amendment.

Common-round phase correlation alone is too limited as a conventional comparator:
it cannot represent different color offsets. Add codebook-blind per-frame
population-template alignment. Starting from common-round alignment, for ten
fixed steps average all 28 aligned preprocessed frames, estimate each original
frame's translation to this template with phase correlation (upsampling 50),
and center the resulting shifts over frames. Use original frames at every step,
so repeated image resampling does not accumulate. Record iterations and runtime.

Give the author model both initializations: common-round and population-template.
Run the original 50-iteration optimizer for each and select by its final retained
image-model loss, never by annotations or imposed shifts. Count both fits in
runtime, save both losses and failures. A failed initialization is recorded;
the remaining successful fit may be used, with that exception disclosed. If both
fail the policy fails. This broadens the strong alternatives before evaluating
whether any of them wins. It is not a new algorithm claim.

Five policies replace four; all other conditions, thresholds, scoring and gates
are unchanged. Save this amendment's hash with every arm and the setting freeze.
The synthetic sign/interpolation control must verify the model's input-coordinate
convention before the real perturbed study. Frozen settings apply to every crop.
