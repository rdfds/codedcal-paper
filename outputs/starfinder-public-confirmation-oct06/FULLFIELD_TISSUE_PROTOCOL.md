# Full-field geometry repair — exposed tissue development

The original frozen tissue adapter fails all three interfaces: 18/47/91 coding
objects and minimum coverage zero. Its per-ROI global registration estimates
up to 73 pixels and 14 stack planes, followed by local flows up to 24 voxels.
The original 96-arm comparison was not run and is not silently completed here.

The author workflow estimates merged-round global registration on the full
field before cropping. This repair restores that spatial information, without
changing candidate domains, decoder cuts, quality gates, books or search core.
It is post-hoc development on exposed tile_1, not independent confirmation.

Use the exact retained full-field photometry LUTs. Form each round's full-volume
maximum over four RNA channels. Estimate unnormalized full-volume DFT alignment
using every Z plane and a fixed factor-two X/Y subsample; refine the integer
translation over its 27 neighboring full-grid shifts using source-image dot
products on a fixed factor-two reference grid. This is a declared resource
port, not exact full-resolution MATLAB DFT. No book labels, imposed offsets or
output scores enter source registration. Round one defines the reference gauge.
Extract the same five ROIs from translated source coordinates, with explicit
axial support; then apply the previously qualified 3D local ROI port, without a
second ROI-level global fit. Record full-field maps and all local flows.

Before real repair outcomes, qualify zero/odd/even translations, support,
coarse/fine signs and the unchanged complete-book actual decoder controls.
Retain every failed original adapter result. Native eligibility is unchanged:
50 coding objects per family over all windows, 10 in designated crops zero/four,
minimum pre-exclusion coverage .90. Failure remains failure; do not lower gates.
A passing repair permits development robustness tests, but cannot retroactively
establish a fresh frozen preparation result or RNA accuracy.
