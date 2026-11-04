# Frozen public tissue confirmation — conditional catalog recovery

The positive repair remains exposed development. STARfinder HeLa culture shows
that a strong multistart image template can tie it; the stronger 2D multistart
still misses BARseq3 color outputs (F1 .791/.611 development, .015/0 spatial
holdout). These results motivate testing conditional robustness, not a promise
of superiority in every preparation. No tissue intensity pixels are exposed.

Source: metadata-selected mouse-brain tissue tile_1 from Zenodo 11176779.
Complete64physicalwords, four rounds, ch00/ch02/ch01/ch03, rotate -90,30planes.
No comparable same-parity absent words exist. Negative-code proxies and RNA FDR
are unavailable, rather than zero. This is another source experiment; all five
windows come from one field and do not count as five preparations.

Declare a resource-bounded Python ROI adapter: full-field per-channel minmax and
64-bin nearest-CDF histogram matching against round1/color1, source -90 rotation;
five fixed128-square ROIs centered at (768,768),(768,2304),(1536,1536),(2304,768),
(2304,2304),all30planes. Per-ROI source-style merged-round 3D DFT followed by the
qualified four-level10-iteration Demons port. This is not full-field MATLAB
pipeline equivalence. Keep source hashes, counts, LUTs, maps and flow hashes.

Conditional target: one attributable code per fully supported source-defined
first-round peak (26-neighbor plateaus >.2 full-field-conditioned maximum,
5x5x3 traces, source per-round L2 norm). Same culture-learned ambiguity cuts:
source .019306540489196777, correlation0, regression0; source exactword lookup,
actual pinned singleshot noise floor .05/cut .775, actual pinned nonnegative
100iterations/PSF0/L1=0/gain fit/cut .2densitymax. One best assignment per object;
uncertain objects abstain. These cuts are not recalibrated using tissue pixels.

Eligibility: all three decoders >=50coding objects across five fixedwindows;
>=10coding objects in both designated calibrationcrop0 and scoringcrop4; each
window >=.90 support before exclusions (empty windows have zero coverage).
This is eligibility of a narrower supported-object target, not passage of the
old .99 whole-interface gates. Report all unsupported/ambiguous objects. Failure
stops the frozen comparison and remains a negative qualification outcome.

Qualify complete four-round/64word loading, synthetic60genes all three decoder
families F1>=.99/no erroneous codes, geometric exact/control inverses, ambiguity
abstention, count-preservinglabel corruption, and all comparator warp conventions
before native pixels. Keep current global_group core byte-identical.

Fit on crop0, evaluate fixed crop4. Initial 2D plane selected only by first-round
supported source peaks. Cases: zero, balanced roundcommon +/-(4,2) with two rounds
each, balanced colors0/1+versus2/3- with identical completeoffsetmultisets; one
out-of-class diagnostic with color shifts (3,1),(-2,4),(1,-3),(-2,-2).
The latter cannot support universal safety of a balanced-group model.

Policies: none; ten-step common; three-start2Dtemplate (zero/round/firstframe,
image residual selection); pinned author barcode50iterations from zero/common/
2Dtemplate (fitcodingloss selection); three-start3Dtemplate with same bounded
image-residual selection; mergedround3Dglobal+local source port applied directly
on scoreobservations; unchanged global_group; exact inverseoracle. All policies
share native preparation, books/cuts; no fitted policy sees imposedmaps or
referencecalllabels. The source port has the extra privilege of score pixels.

Before outcomes freeze observable rule: among conventional policies choose the
lowest coding loss on a first-round-selected scoringplane; choose global_group
only if its score loss beats the conventional minimum by >=.005. Neither known
offsets nor output F1 enter the decision. Rules are evaluated, never assumed safe.

96arms = four cases x eightpolicies x threedecoders. Primary identityF1 at voxel
radius3, sensitivity1.5/5, maximum-cardinality/minimum-distance matching; retain
compositionTV, yields, abstentions, support, undefined-emptyrefs, all failures and
per-policy runtime. Oracle reuse only after exactcandidate+kerneltrace equality.
Reference is native algorithm output, not independently known RNA truth.

Confirmatory interpretation: primary zero/common/color rule avoids >.02 F1 loss
relative to best conventional policy, and finds >=.05 gain when one exists.
A tie supports conditional robustness but not distinct method benefit. Outclass
harm blocks general recommendation claims. No Nature Methods confidence upgrade
from native eligibility alone; all previous failures/conditional temporal limits
remain. This specification/executable/control hashes must be committed before
pixel preparation. No post-outcome repair is confirmatory on this field.
