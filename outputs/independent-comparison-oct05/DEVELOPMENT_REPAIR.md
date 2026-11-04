# Development-only registration repair before validation

The initial population-template baseline failed on unperturbed development data:
its mean-centered shifts had RMS 45.89 pixels and maximum component 84.91 pixels.
Correlation-decoder identity F1 fell from .8643 to zero. A BARseq readout contains
only a subset of the sources; unit-normalized phase correlation across different
subsets amplified weak/noisy frequency components and admitted a distant peak.
The synthetic example did not expose that preparation-specific failure.

The stopped run and all completed/partial arrays are retained in
`development-unbounded/` and `development-unbounded.log`. No validation crop was
opened. No detection threshold, scoring endpoint, dose or qualification target
is changed. This is an adapter/baseline repair, not a favorable-outcome selection.

For the ten-step codebook-blind template baseline, replace unbounded phase
correlation by zero-mean amplitude linear correlation via FFT, search lags only
within [-10,10] pixels on each axis, and use local quadratic subpixel peaks.
The ten-pixel bound matches the author model's preset search padding and is
fixed for every image, not chosen from the known imposed shift. The actual
model does not strictly constrain its returned shifts by its padding parameter;
record its entire output. Common-round phase correlation remains unchanged.

Both model initializations remain, selected by model loss. The observable
preflight now uses the repaired template's corrections, still independent of
annotations and imposed shifts. Re-run all development conditions and the
known-source sign/recovery control. The decision threshold must freeze on that
complete development run before validation. Failures remain outcomes.
