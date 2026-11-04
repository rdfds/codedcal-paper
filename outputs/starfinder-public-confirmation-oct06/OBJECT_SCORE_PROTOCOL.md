# Direct input-object score — exposed diagnosis

The tissue pixel-loss rule fails: additional nonrigid correction lowers the
warped-image coding loss while losing 6–10 F1 points; the color search restores
calls but its pixel-loss advantage is below .005. Diagnose a target/representation
mismatch before changing thresholds. This is development, not confirmation.

Track first-round source peaks by their observed color, using fixed x/y32:96,
z6:24, threshold .2 and complete support. Carry each anchor by its own fitted
first-round map. Evaluate each workflow by sampling the unchanged input voxels
at its inverse map, with integer nearest sampling and unique voxel accounting
for nonlinear maps. Do not evaluate newly interpolated intensity images. The
shared inputs are already source-conditioned/deconvolved; do not call them raw
camera photons. Penalize unsupported objects rather than dropping them.

Use 5x5x3 sums, per-round L2 normalization and the maximum projection on the
normalized complete physical book. Each anchor receives equal weight; missing
readouts receive maximal loss. Labels/reference calls/imposed maps are excluded.
Nonrigid flows are replayed only to reconstruct inverse sampling coordinates,
and their hashes must match the original retained diagnostics.

Abstain with fewer than ten tracked anchors. Keep the .005 margin. Include no correction as a conventional workflow. Among
conventional scores within .005 of the minimum choose the simpler policy in
fixed order: none, common,2Dtemplate,3Dtemplate,author,source-local. Choose group
search only if it beats the conventional minimum by >=.005. Report bootstrap
uncertainty as a diagnostic, not an RNA correctness guarantee. Retain every
case/failure/abstention, including out-of-class harm. Numerical controls must
cover per-round gain invariance, exact correction, unchanged input sampling,
missing support and source-coordinate conventions before native score outcomes.
A successful development rule still requires a new frozen field/acquisition test.
