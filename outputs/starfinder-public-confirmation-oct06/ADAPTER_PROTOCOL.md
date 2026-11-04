# STARfinder source-family adapter — development qualification, 6 October 2026

Before opening intensity pixels, use only metadata-selected culture Position351.
Other culture fields and tissue intensities remain unopened. The complete source
book has 998 coding words and 26 same-end-base-parity absent words. Decode all
1024 parity-valid words, calibrate with 998 physical words, and scale negative
counts by 998/26. These negatives are an operational proxy, not RNA FDR.

Read all thirty planes of each uint8 RNA stack, rotate by metadata -90 degrees
(clockwise, no resampling), retain source color order ch00/ch02/ch01/ch03.
Compute per-channel min/max and histogram profiles using every source voxel.
Normalize to uint8 and histogram-match to round-one/color-one with a declared
64-bin nearest-CDF LUT port. MATLAB's implementation/tie conventions are not
certified identical. No source intensities or ROI outcomes select a LUT.

Register full 3D merged-channel maxima with unnormalized Fourier correlation,
as DFTRegister3D does. Integer shifts and zero-filled application are checked
against a direct reference. Array axes are z,y,x; source shifts are row,col,z.
Source default rotation and complete rounds are retained. No maximum projection,
downsampling or imposed truth enters registration. Use scipy real FFTs and one
channel at a time to fit the local 8 GB machine. Full source files remain.

This initial variant uses global registration. Optional source-local MATLAB
Demons is not implemented or claimed equivalent. If a later local variant is
needed it must have separate controls/freeze and be called development repair.
That variant and source local registration must be included among strong
alternatives before an independent comparative conclusion is earned.

After full-volume photometry/geometry, extract five fixed 128x128 windows at
(y,x)=(374,374),(374,1122),(748,748),(1122,374),(1122,1122), retaining full depth.
No content-selected projection/plane/ROI is allowed. Decoder data use z=13..15
(zero-based), y,x=16..111; score central plane z=14, y,x=24..103. The domain is
explicitly a central-plane endpoint in a 3D-conditioned acquisition, not full
volume performance. Keep support through all global shifts.

Three decoder interfaces: (1) source-family 3D regional maxima on each reference
 RNA color, adaptive intensity cut .2 of the full reference-channel maximum
 (1 after the declared LUT), not a cropped maximum; connected-plateau centroid rounded to a
voxel, 5x5x3 summed intensities, unique dominant color per round and exact word
lookup; (2) pinned BarDensr singleshot correlation (.05 noise, .775 peak cut);
(3) pinned nonnegative regression, 100 iterations, zero L1/PSF, codebook gains,
.20 of density maximum. The source q-score parameter is declared but never
applied by FilterReads; do not invent its use here. Ports and kernels are not
claimed as exact full MATLAB pipeline execution.

Numerical qualification: full source semantics; histogram monotonicity and
independent LUT reference; rotation/axis/FFT known-shift controls; 3D regional
max/voxel-color calling; native kernel recovery of 60 synthetic genes drawn
before pixels from the complete source book, with zero absent calls; support
and exact integer inverse controls. Freeze source/core/runtime hashes after
controls, commit adapter before any native execution.

Native gate per decoder: all five windows completed; >=50 coding calls, scaled
absent-word proxy <=.10, central scoring support >=.99. Report every decoder's
status separately; all three must pass for this declared three-family interface
to qualify. No threshold relaxation or perturbation comparison after failure.
Further development may repair the adapter but requires a fresh reserved test.
This specification is native qualification only. It does not freeze the full
comparative study, create independent RNA truth or validate the temporal target.
