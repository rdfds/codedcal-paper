# Fresh culture adapter qualification, before comparative outcomes

The fixed central 256-pixel, middle-plane patch is the sole fresh input. Before
RNA decoding, register each published frame6 fiducial to round1 by bounded image
correlation (32-pixel XY search, quadratic subpixel peak), canonicalized to round1.
Apply that same measured correction to both RNA colors of its round. This is a
native source-adapter operation using an independent image channel, not injected
challenge truth or a scored RNA catalog. Retain maps, support and post-correction
fiducial residuals. No data-based crop, threshold or method selection is allowed.

Then use the pinned sigma3 background subtraction, nonnegative clipping, sigma.7
smoothing, and per-readout positive 99.9-percentile scaling. The decoder context
is the fixed central 224-pixel square (16:240), score/coverage 24:232. Use the
same actual BarDensr correlation (.775, noise .05) and 100-iteration regression
(.2*max density) kernels. This is an explicitly qualified adapter, not a claim of
reproducing the complete Serval or MERlin author production pipeline.

Native/no-effect/integer-inverse/fractional operator gates must pass, with >=50
nonblank calls and >=99% scoring-domain support. Also retain native blank proxy
and its one-sided Poisson upper95 bound. For fresh-source confirmation eligibility,
require upper95 <=.1 if the source codebook provides named blanks; no arbitrary
normalization of the proxy is permitted. Failing this gate is a source failure,
not a successful independent confirmation. Qualification may repair demonstrable
mapping/implementation defects before comparative freeze; never tune against
comparative outcomes or silently change the frozen thresholds.
