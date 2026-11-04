# Observable preflight decision — frozen before nonzero outcomes

Only three unperturbed development arms had been printed when this extension was
written. No nonzero-condition outcome and no validation outcome has been read.
The original comparison proceeds regardless of whether this extension qualifies.

The practical decision is common-round versus per-frame population-template
alignment. The observable feature is RMS of the difference between their 28
estimated corrections. It uses only the observed images and can be measured
before decoding or using annotations; it does not use the imposed displacement.
This is a preflight that computes both registrations, not a zero-cost diagnostic.

A positive target is per-frame alignment improving primary identity F1 by at
least .05 over common-round alignment, for each decoder separately. Train a
one-threshold rule on the 12 nonzero development conditions, from the fixed grid
[0,.25,.5,.75,1,1.5,2,3,4,8] pixels. Select maximum balanced classification
accuracy, then minimum false positives, then higher threshold. Require at least
three positive and three negative conditions and balanced accuracy >=.80 to
qualify. Save the selected threshold and all development outcomes before loading
the validation crop. Do not fit a regression or add features after seeing that
crop. If qualification fails, retain the failure and make no predictive claim.

For validation report confusion counts, balanced accuracy and the selected
workflow's identity F1 alongside always-common, always-per-frame and the author
barcode model. A false-positive call is counted even if it is cheap. Report cases
where selected alignment harms identity by >=.05, and runtime including the
preflight. Report decoder-specific rules; two perturbation seeds are technical
replicates. This does not establish independent acquisition/animal transfer
between the two BARseq crops. It supplies an independently annotated public
assay, alongside the existing MERFISH consistency evidence, with distinct scopes.

Neither a successful rule nor a successful comparative result certifies identity
from images alone: codebook-label permutation remains an unobservable ambiguity
without an independent reference. This preflight is empirical within its tested
scope, not a general identity guarantee.
