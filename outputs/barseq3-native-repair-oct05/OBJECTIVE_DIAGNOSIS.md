# Registration objective diagnosis (6 October)

This is a post-hoc diagnostic on the two exposed BARseq3 fields, not independent
validation. The primary recovery study and interior-fit sensitivity remain frozen.
No reserved field, private data, or paid job is used.

Before examining diagnostic outcomes, compare all six field/case combinations
(fields 0/1; zero/common/color) using the pinned author single-code objective:

1. Evaluate no correction, existing full-fit and interior-fit model corrections,
   and known inverse injected correction. Record the original padded objective
   and a common fully observed 80-pixel interior. The latter uses the same
   interpolation/code normalization but fixes denominator and spatial support.
2. Start the author's unchanged 50-step optimizer at the known inverse injected
   correction on the full 128-pixel input. This truth-assisted run is a mechanism
   control only. Decode its rounded correction with both existing frozen kernels,
   windows, cuts and three matching radii. Record native-call consistency,
   composition, blanks, yield, support and displacement from the oracle.
3. For each color case, evaluate nine equally spaced geometries between the
   existing failed/full fitted map and oracle. This is an objective profile, not
   a search or a deployment initializer. Retain gradients at the oracle.
4. Independently verify single-code energy using NumPy and integer shifts; check
   that the oracle exactly restores the fixed interior. Retain every failure.

Interpretation: a lower oracle loss with successful oracle-start optimization
supports an initialization/local-basin limitation. Preference for the failed map
or harmful oracle-start drift supports an objective limitation. Both may occur.
Small fractional drift alone does not establish a harmful objective. Native
conditioning, crop context, photometry and overlap remain possible contributors;
a single-code assumption alone does not demonstrate an overlap mechanism.

The known inverse is not available to a user and cannot qualify a recommendation.
A repair will be chosen only after this diagnostic; its design will be recorded
separately before scoring. All improvements on these fields are development.
