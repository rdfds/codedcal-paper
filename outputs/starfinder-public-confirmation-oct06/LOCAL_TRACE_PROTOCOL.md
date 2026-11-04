# Source-family local 3D registration and normalized traces — development

Retain all earlier gates. Full-depth global-only qualification completed 15/15:
source 195 coding/1 negative (scaled proxy .197); correlation 17/0 (yield miss);
regression 1440/26 (proxy .693). None passes the unchanged native criterion.

The source workflow runs local 3D Demons after global registration (10 iterations,
AFS 1, floor(log2(30))=4 pyramid levels). This port uses SimpleITK symmetric voxel
geometry, levels 8/4/2/1, 10 iterations each, field smoothing 1. Unlike the MATLAB
full-field implementation, local fields are fitted on the five existing 128x128
full-depth windows. This is a declared ROI/runtime port, not exact author-pipeline
equivalence. Preserve flow hashes, interpolation/support, skipped identity cases
and all native outcomes. Never use imposed truth to fit source-local registration.

Keep full-depth source-defined first-round candidates and 5x5x3 extraction. The
author's ExtractFromLocation divides each round's four-vector by its L2 norm
plus 1e-6 before assigning its dominant color. Use that normalized four-vector
for all trace interfaces rather than introducing a mean-intensity scaling with
a noise floor borrowed from pixels. Source exact-word assignment is unchanged;
actual correlation (.05/.775) and nonnegative (100 iterations/.20 maximum)
parameters and full 998+26 book remain. This targets relative identity evidence,
not calibrated molecular abundance. Multiple numerical gene hypotheses remain
explicit; no confidence threshold or blank ceiling is tuned to win.

Native eligibility still requires all five windows per family completed, >=50
coding calls, scaled same-parity-negative proxy <=.10 and support >=.99. Support
is the fraction of otherwise in-domain first-round candidates with complete
5x5x3 support in all frames; unsupported candidates are reported as abstentions,
not hidden or assigned identities. This object-domain definition differs from
the previous whole-volume pixel support and is explicit. Empty windows give no
calls. All three interfaces must qualify before comparative testing.

Numerical controls precede this variant's outcomes: identity local flow;
pull-warp x/y/z sign and support; known integer/3D global inverse; direct
per-round normalization against author arithmetic; fresh native kernels with
full-depth synthetic genes and zero negatives. Position351 is exposed development.
Reserved fields/preparations remain unopened; the independent study is not frozen.
