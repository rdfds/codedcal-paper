# One-field native qualification retrieval, pre-comparison

Only the nine public 4T1-culture TIFFs named `merFISH_merged_01_000.tiff` through
`merFISH_merged_09_000.tiff` in the retained released Harvard record are eligible.
No other FOV or preparation is eligible. All nine original whole-file MD5 values
and sizes must be pinned from that metadata before network I/O. Each TIFF is
232,717,320 bytes; total transport <=2,094,455,880 bytes (nine requests), 180-second
per-file timeout, ten-second spacing, one reader, no automatic retries. Stop on
HTTP429/Retry-After, length/hash mismatch or disk reserve <1 GiB. Keep only one
whole temporary TIFF at a time. Do not erase partial failures or receipts.

Check all publisher MD5s and record SHA256/response headers before extracting
fixed pixels. Readout frame numbers use the publisher organization, zero-based
indices, with microscope flips/transpose recorded and consistently applied.
Qualification input is the fixed central 256x256 patch, z index 7 (the middle of
15 documented planes), all 16 RNA readouts from rounds1..8. Round9 is retained as
native metadata/context but is not an RNA readout in the 16-bit codebook. Extract
also each round's frame6 fiducial for the same patch, to document rather than
silently substitute fiducial alignment. Native preprocessing is fixed before
outcomes: per-readout Gaussian background subtraction sigma3 and nonnegative
clipping, Gaussian smoothing sigma.7, then each frame divided by its positive
99.9th percentile. No fitting of comparative calibrators or challenge outcomes
is permitted in this reader. No crop/count/threshold search. Nonqualifying native
counts/support or invalid mapping means source qualification fails.

Only after checksum verification may the redundant original temporary copy be
removed. Keep extracted arrays and per-file extraction parameters/hash receipts.
The released record is one preparation; patches/fields would not create extra
independent preparations. Access and native qualification alone are not an
independently confirmed comparative benefit.
