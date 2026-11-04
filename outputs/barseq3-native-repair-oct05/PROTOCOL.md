# Native-conditioning repair, before full-field image processing

The previous partial-adapter external gate failed and is closed in
../barseq3-transfer-oct05. Do not overwrite it or call its rule externally tested.
Both central native crops are exposed. This repair uses the published processing
recipe rather than tuning a detection cut until the existing gate passes.

Source: Zenodo 19024678, adult_cb_raw_FOV20_FOV40_11rds (1), 934,285,696 bytes,
MD5 a03d84ad28824b9fef6ff2ffd514eee3. Author code commit
32df2c4793b8ab0534b73cc56a506ae16d4048a8 and supplied gene-only configuration.
Download only this one bounded file, with exact size and full MD5 verification;
the earlier 128-MiB transport cap applies to the closed selective attempt.
The new cap is the source's exact size plus metadata. Require >=3 GiB free
before download, process only one dense job at a time, and retain raw input.
No cell-expression/annotation values or other huge archive is needed.

Native recipe: first ten rounds and G,T,C,A. Fit the authors' camera homography
from channel 3 to 2 per round using their blob/SIFT features, 90% descriptor
matches, <10-pixel point distance and 1-pixel RANSAC. Apply it to channels 1/3.
Fit the authors' partial affine round-to-round maps from all-channel summed
images to round zero, retaining 80% of descriptor matches with blob/SIFT and
1-pixel RANSAC. Apply the published
spectral inverse. Use the authors' common-overlap bounding rectangle and 20-pixel
margin. No translation fallback, bad-point trimming, adaptive crop or changed
RANSAC criterion will rescue a failed native geometry fit.

Two-stage percentile normalization uses precision 0.005 and sigma-20 background
subtraction, followed by the source Demons recipe (150 iterations, smoothing
2.0, round-zero reference). The full release does not provide the configuration's
five normalization FOVs, so source-zero alone supplies the fixed calibration
percentiles for both fields. Disclose this restriction; it is not the authors'
complete multi-FOV/cell-segmentation production run. Stream per-round operations
and restrict only the final native resampler output to the ROI, rather than
keeping multiple gigabyte stacks. Gaussian and crop-restricted resampling must
be numerically checked against the native operations before use. Source audit
before image processing corrects two earlier assumptions: the overlap mask
requires nonzero unmixed signal in every frame/channel, not round sums, and
the pinned BarDensr background routine overwrites intermediate blurs, so only
the last-axis sigma is applied. Preserve this implementation behavior rather
than silently substituting a two-dimensional filter. The upstream dependency
pin is dfdafe9055be322b286369d21e3349204c1f6fc3; the release does not pin its
installed dependency version, so this remains an explicit dependency assumption.

Default correlation uses author cut 0.775 and noise floor 0.05. The separately
specified actual regression kernel retains 100 iterations, zero PSF/L1, estimated
frame gains and cut 0.20 times maximum. Preserve every native call/control and
shared-sequence ambiguity. These settings are fixed before new decode outcomes;
the failed 0.75/0.01 partial adapter is not retrospectively replaced.

Development ROI: source zero, full-field [1088,1216) in both axes. A spatial
validation ROI, source one [1408,1536), is fixed now, outside the earlier raw
context. Its full-field geometry may use unlabeled images, as required by the
workflow. Do not decode it before source-zero qualification and a committed
normalization profile. Native gates remain >=50 nonblank interior calls,
adjusted blank proxy <=0.10, finite densities and >=99% all-frame geometric
support. Interior is the central [24,104) of each ROI. Preserve failures;
no third crop is chosen after outcomes.

This is adapter development on an exposed external preparation, with a later
spatial holdout, not untouched-preparation confirmation or RNA truth. Successful
native qualification only permits a separately frozen recovery comparison; it
does not itself establish calibration benefit. The CSHL two-pixel decision cut
is not refit. Any controlled recovery must distinguish synthetic offset truth
from native-call consistency, include the stronger model-based alternative,
and retain all failures. Neither study validates a physical randomized readout
or changes the failed ISS/S3 and two estimator-modality-class qualifications.
