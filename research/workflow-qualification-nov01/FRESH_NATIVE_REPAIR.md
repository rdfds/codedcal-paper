# One bounded source-conformant native-adapter repair

The first native attempt is preserved:3 correlation/12 regression nonblank calls,
with inadequate blank-proxy precision. It used an improvised99.9-percentile
intensity scale and omitted the source's deconvolution. Source inspection shows
Serval's actual standard pipeline uses high-pass sigma2, native10-iteration
Richardson-Lucy (window9,sigma1.2), low-pass sigma1, and histogram90th-percentile
initial scales. The magnitude/noise gates assume those intensity units. This is
a demonstrable adapter mismatch, not evidence that a different calibration wins.

One fixed repair is permitted before any comparative fresh outcomes: use those
exact existing native filtering/scaling functions and the actual source's
high-pass sigma3/upsample100 fiducial alignment. Keep the SAME middle plane,
256-pixel patch,16-bit literal codebook, call thresholds,50-call minimum,99%
voxel-support gate and blank-proxy upper95<=.1 eligibility rule. Preserve all
old results. No crop/threshold/parameter grid or further method choice is allowed.
Native controls only, both actual BarDensr backends; no comparative policy/challenge
outcomes may be opened. This is a qualified source-conditioning adapter, not a
claim to reproduce the complete published Serval decoding/optimization pipeline.
If still unqualified, retain that failure; do not force a confirmation claim.
