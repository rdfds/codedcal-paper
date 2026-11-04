# Interior-fit sensitivity

This additional diagnostic is declared after source-zero zero/common recovery
outcomes and the start of its color case. It is post hoc, not another independent
confirmation and not a replacement for the frozen 60-arm comparison. Neither
that running implementation nor its crops, policies, cut, offsets or outcomes
are changed.

Potential confound: synthetic shifting of a 128-pixel ROI fills outer pixels
with zeros. Those artificially displaced edges may help a calibration fit.
Test every original field/condition/policy/decoder again with fitting restricted
to [24,104) interior pixels. All imposed shifts leave those input pixels real.
Ten-step templates are formed from warped full inputs, but correlations use only
the interior, with explicit support checks for every template iteration. Phase
initialization uses the same interior; the actual author model fits that interior
from its same two image-only initializations, 50 steps. The calibration estimator
and its comparison domain, not its declared optimizer/cut, change. The smaller
fit domain is a declared sensitivity that can itself reduce information.

Apply the resulting rounded maps to the unchanged 128-pixel observed arrays;
use the original common 96-pixel decoder window, native references, gene-class
endpoint, three radii, point-proxy/support gates and >=2-pixel rule. No tuning,
other crop or favorable subset is allowed. Two fields × three cases × five
policies × two kernels = 60 arms. Count all failures and missing arms.

If a corrected decoder input window exactly equals a retained primary-policy
window, reuse that policy's calls only after bytewise array equality. Otherwise
execute the actual kernel. Retain the reused input/source-result hashes and
actual execution cost separately; this is not an independent kernel repeat or
a faster scientific workflow. Oracle windows must still reproduce native calls.

Require the original comparison to complete before this sequential job starts;
use one numerical process at a time. A blocked primary gate blocks the follow-up.
Report both cohorts. A disagreement identifies a limitation; it is not permission
to discard the first cohort or soften a threshold. Native-output agreement is
not independent RNA truth, and marginal RMS matching does not measure temporal T.
