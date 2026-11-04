# Retained-observation comparison of arithmetic and geometric temporal contrasts

17 October 2026. This is a separately specified, retrospective extension of
Source Data 70, whose outcomes are already exposed. Freeze this protocol and
the implementation before calculating new comparative outcomes. No new source,
simulation realization, reserved pixels or physical intervention is introduced.
Original Source Data 69/70 files, gates and failures remain unchanged.

## Complete cohort and targets

Use all 126 completed Source Data 70 worlds: nine views of eight acquisitions,
both noise laws, four primary motion laws and three violations. Retain every
field, both independent code assignments sharing a physical field, the
known-source comparator, observed group counts and the 512-prefix sensitivity.
Verify each retained array against its existing state-file SHA-256 before use.
Stream inputs; do not copy observations or generate new realizations.

For per-group window contributions (C_A,R_A,C_B,R_B), compare:

- Arithmetic: (E C_A + E C_B)/(E R_A + E R_B), the shipped primary contrast.
- Geometric: sqrt((E C_A/E R_A)*(E C_B/E R_B)), the existing window-split
  construction. A physical positive-attenuation call requires both window
  ratios and both observed reference means positive. Two negative window
  ratios must not earn a positive call merely because their product is positive.

For each, calculate its own physical baseline and actual-diagonal observable
target from the retained population moment table. The geometric baseline
target is the square root of the product of the two baseline window ratios;
the arithmetic target is the original pooled baseline ratio. Report their
difference explicitly. They coincide for homogeneous/common T but need not
coincide under source-dependent motion with different window weights. A
change of target is not improved recovery of the original arithmetic target.

For a source-uniform multiplicative window-B kernel factor f>0, window ratios
become f*a and a/f under the common-motion model. Their geometric mean is a;
the pooled arithmetic ratio generally changes. Verify this algebra in controls.
Source-dependent changes and code-dependent geometry have no such guarantee.
Actual-diagonal targets do not represent actual randomized-code expectations
when code exogeneity fails; retain the resulting discrepancy rather than
claiming the diagonal is independently known physical truth.

## Confidence sets fixed before outcomes

Arithmetic uses the shipped pointwise 95% Fieller set for paired group means,
with the same usability rule as Source Data 70. Recount arithmetic set types,
point estimates and coverage against its archived field records.

Geometric uses two pointwise 97.5% Fieller sets, one per window. Intersect each
full set with positive ratios, multiply the positive sets and take the square
root, retaining every interval component and unbounded endpoint. Bonferroni
gives nominal joint coverage at least 95% if both marginal approximations are
calibrated; independence between windows is not required. Calibration remains
approximate and can fail for concentrated group contributions. A usable call
requires both original window sets bounded and wholly positive, with positive
observed reference means. A finite transformed subset cannot rescue an
unsupported window. This conservative rectangle may sacrifice precision.

No delta-method interval is substituted after seeing outcomes, and no
data-selected choice of variant, frequency, confidence rule or source subset
is made. This compares the existing point constructions with declared
uncertainty treatments, not an optimality theorem for every possible interval.

## Endpoints and qualification

Separately per variant, method, scope, view, noise law, world and frequency,
retain point availability; positive/negative/nonfinite support failures; full
set types; numerical usability; physical-target and actual-diagonal coverage;
usable-and-covering frequency; finite-ratio sample bias/RMSE/tails; interval
width; and usable coverage with upper bound <1 (small positive-motion resolution).
Report the geometric error relative to the original arithmetic baseline too,
without qualifying it against an incompatible target.

Apply the existing descriptive gates to each variant's *own* physical target:
absolute finite available sample-mean bias <=0.01, >=90% numerical usability,
>=90% joint usable coverage, and a 2,000-draw whole-physical-field bootstrap
interval for joint usable coverage including nominal 95%. Preserve paired code
assignments. Retain above/below-nominal failure direction; all-success bootstrap
degeneracy is not undercoverage. These are inherited study criteria, not new
field standards or proof of exact coverage. More passes alone do not establish
a better estimator: compare width, resolution and changed targets alongside.

Retain the original D computation without adapting its cutoff. For every
variant, count available absolute physical-target errors >0.01 with quiet D,
and invalid diagnostics. Quiet D never certifies either variant. Report all
three violations, not just the source-uniform case favorable to geometry.

## Controls, execution and reporting

Before outcome calculation, test exact uniform-factor cancellation, a
source-dependent counterexample, heterogeneous-weight target differences,
all confidence-set topologies, product-set endpoint propagation, strict JSON
and geometric point rejection for negative/unsupported windows. Check arithmetic
agreement with the frozen function on retained examples. Freeze protocol,
implementation, source manifests and controls in Git before outcome execution.

One bounded local worker streams existing arrays and writes compressed per-case
records, complete summary tables, source hashes and a resumable per-world ledger.
No replacement scientific jobs, new simulations, paid compute or outcome-driven
repair. Implementation errors require a separately recorded repair; do not
silently replace a failed attempt. Independent scalar/set-membership recounts
must pass before manuscript integration. Preserve negative results and update
the main figures only to the extent this conditional comparison supports them.
