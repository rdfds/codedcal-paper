# Independent public BARseq comparison — frozen before image outcomes

## Question and scope

Does a strong demixing decoder remove the need for readout calibration, and
does common-round alignment suffice when color channels have different offsets?
Separate recovery of positions from recovery of barcode labels. A useful result
must concern the required workflow, not claim that our estimator is the winner.

The public BarDensr benchmark is a preparation independent of our MERFISH
archives. Its two annotated crops are not two independent preparations or two
known independent animals. The authors' annotations are a hand-curated reference,
not physical molecular truth. Added shifts supply perturbation truth only. This
study cannot validate randomized readouts or temporal variance separation.

## Exposure and sources

Before this protocol, only repository source, notebook code (not outputs), HDF5
schema/attributes and file hashes were read. No benchmark image or annotation
values were decoded or scored by this project. The source checkout was already
present for unsuccessful MERFISH qualification experiments; those failures remain.

Source: <https://github.com/jacksonloper/bardensr>, commit
`dfdafe9055be322b286369d21e3349204c1f6fc3`.
Author paper: <https://doi.org/10.1371/journal.pcbi.1008256>.

Development: `ab701a5a-2dc3-11eb-9890-0242ac110002.hdf5`, SHA256
`d78bd29c10642c7f0ee30ae2d29445e0d0b5c44a6329846fc926f0f90751fbf6`.
Held-out crop: `bc_L2_MAX_M1-L2_000_002_x1500_y600_v004.hdf5`, SHA256
`520cabe0d7ac2cb7ec869219a7a477c846f69faea3d59c38090ab6048b793b72`.
The registration and colorbleed HDF5 variants are not additional preparations.

## Frozen experiment

Use every good author annotation and every predicted call inside the common
interior x,y in [24,176). Report excluded annotation statuses and border counts.
No adjustment of the annotations, registration gauge or scoring region from
outcomes. All registration uses images and the released codebook, never labels.

Conditions: unperturbed; common-round translation; common-round translation plus
static color offsets. Nonzero RMS doses 1, 2, 4 pixels, seeds 501 and 502. Draw
normal translations and center separately over rounds/colors, then scale each
component to its specified 2D RMS. Common-round dose applies to all four colors
in each of seven rounds. In the color-offset condition, add independent static
color offsets of the same RMS. Actual total RMS is recorded, not equated between
these two different perturbation structures. No image corruption is real drift.

Four calibration policies, all counted: none; common-round phase correlation of
the sum of four preprocessed channels (50-fold subpixel upsampling, round zero
reference, mean-centered gauge); author's barcode-conditioned model registration
(maximum wiggle 10, 50 fixed optimization iterations, initialized from common
round alignment); oracle inverse added translations. The last is an interpolation
control, not an available workflow. Model failure is reported, not silently
replaced. Apply corrections with scipy linear shift, constant zero outside input;
all conditions share the interior region. Record correction RMS against known
inverse perturbations, but never use that truth to choose a correction.

Two distinct decoder families: the author's correlation evidence kernel
(`estimate_density_singleshot`, noise floor .01) and nonnegative regression
(`estimate_density_iterative`, 100 iterations, zero PSF, zero L1 penalty,
estimated frame gains, other options default). No nonzero-PSF axis repair is
needed or used here. These are actual pinned kernels; this is not a claim of
reproducing the complete MERlin/Serval/BarDensr production pipeline.

Preprocessing follows author tutorial: per-frame minmax, positive Gaussian
background subtraction sigma (0,10,10), then per-frame minmax. For registration
use that preprocessing on observed inputs. For decoding, warp the observed
inputs, then perform the same fixed preprocessing. Fitted decoder gains are part
of each complete workflow and their runtime is counted. They cannot be described
as frozen gains or as a prediction made before observing the image.

Only unperturbed development data may set detection thresholds. For correlation
choose among [.50,.60,.65,.70,.72,.75,.80,.85,.90]; for regression choose among
[.03,.05,.075,.10,.15,.20,.30,.40] times that density tensor's maximum. Choose
highest label-independent-matching identity F1 at radius 5; tie chooses higher
threshold. Peak neighborhood is (0,1,1). Freeze selected settings before any
nonzero condition or the other crop is decoded. Report the full development
threshold curve and half/double thresholds as secondary sensitivity on retained
densities. The threshold is never refit for a policy, dose or validation crop.

## Endpoints and qualification

Primary: barcode identity F1 using maximum-cardinality, minimum-distance
one-to-one spatial matching, independent of labels, at author tutorial radius
5 pixels. Also report spatial F1, wrong-label fraction among spatial matches,
composition total variation, reference/prediction counts, candidate ambiguity,
and radii 1.5 and 3 as fixed sensitivities. Also compute label-constrained
one-to-one F1 to expose matching dependence; it is not substituted for primary.
No blank-derived FDR: the released book has no qualified blank class.

Qualification for opening the validation crop: both unperturbed development
decoder identity F1 >=.70 with >=50 interior annotations, at least 50 predicted
calls, finite density, zero skipped spatial components. If it fails, repair only
on development data, record protocol revision before validation, and retain
failed runs. A numerical/adapter repair must not change the scientific gate
quietly. Recovery must exceed zero perturbation controls, and oracle must not
destroy unperturbed performance by >.05 F1 on either decoder.

Controls before validation: exact-reference score=1; count-preserving cyclic
label corruption degrades identity but preserves spatial F1 and composition;
known 20-pixel geometry translation degrades spatial matching; inverse restores
it. Synthetic overlapping two-code sources distinguish regression from a
single-vector assignment conceptually; this does not qualify a real preparation.

Save all calls, densities (float32 compressed), applied transformations, hashes,
per-stage runtimes and failure tracebacks. Results are paired within each crop;
two seeds and different shifts are technical perturbations, never biological
replicates. No population confidence interval from these two crops. Report every
arm and unsuccessful policy. Stronger workflows winning is a valid outcome.

## Interpretation rule

A reproducible workflow requirement would strengthen a comparative paper. Merely
replicating BarDensr's published advantage under mixing is not a new result. A
simple effect of injected shifts does not establish a deployable decision rule,
independent molecular truth, independent prospective confirmation, or Nature
Methods acceptance. Do not promote manuscript claims until the external result
survives strong calibration and decoder alternatives and these limitations.
