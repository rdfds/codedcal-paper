# Integer-grid sensitivity, specified after completed benchmark outcomes

The original digital perturbations and corrections use linear interpolation.
Some oracle outputs improve on the unperturbed annotation agreement, so those
outputs are not physical accuracy ceilings. This additional diagnostic removes
interpolation blur from both the imposed perturbation and every applied
correction. It is a post-hoc sensitivity on both exposed annotated crops, not
new independent validation or a replacement for the original results.

Fixed conditions per crop: native zero; common-round integer displacement at
amplitudes 2 and 4; common-round plus static-color integer displacement at those
amplitudes. Seed 504. Opposite integer pairs plus one zero round preserve the
global coordinate gauge; two opposite color pairs preserve the color gauge.
Amplitude describes the coordinate range, not exactly the component RMS.

Fit phase, ten-step strong common-round template, ten-step per-frame template
and the actual 50-step author barcode-model optimizer on each observed image.
Use the original absolute-coordinate bounded-template update for both ten-step
fits, rather than changing them to cumulative residual updates. This source
review correction is made before any integer-control fit or outcome is run.
Choose the model's two image-only initializations by its own final objective.
Round every applied fitted correction to the nearest integer pixel. The
frozen decoder thresholds, actual two decoder kernels, three matching radii,
author annotations and image-derived >=2-pixel rule remain unchanged.
Normalize native images once using the original pinned preprocessing, before
any perturbation. Fit and decode those prepared images without another global
normalization. Decode the same [16,184) window for every arm, then restore full
image coordinates before scoring [24,176). This window remains inside exact
oracle support even for the largest added integer offsets; it prevents lost
outer pixels from changing regression gains or the density maximum under the
oracle. The cropped decoder domain is a declared sensitivity, not the primary
200-pixel benchmark. Verify its oracle calls equal its zero/no-correction calls.
Include no correction and inverse integer offsets. Count all six policies,
two families, five conditions, two crops: 120 planned arms, with all failures.

Assert inverse-added offsets reproduce native raw pixels exactly throughout
the scored interior. Report residual geometry and annotation agreement; do
not select a favorable decoder or displacement after results. This stress test
can resolve a resampling confound but cannot establish preparation independence,
a prospective recommendation or molecular truth beyond the author annotations.
