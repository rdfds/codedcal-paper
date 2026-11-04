# Controlled image/codebook workflow challenges

This separate opt-in command preserves the historical estimator CLI byte-for-byte and tests computational catalog preservation under known
digital translations. It does not validate physical motion, RNA identities,
clinical decisions or the prevalence of naturally occurring errors.

`python -m codedcal.workflow_cli input.npz --output plan` writes the exact challenge
plan without importing an optional decoder or executing workflows. `--qualify`
runs native, no-effect and integer-inverse controls and checks the fractional
composed operator. `--run` explicitly runs the qualified challenge comparisons.
The default is the mixed conventional suite: one round, one color, one
round/color interaction and one arbitrary-readout probe. A preflight projection
check inspects the actual maps before importing an adapter. The fourth probe
must add a nonadditive direction beyond the interaction probe. `design.json`
and `design.txt` report representative families, sampled rank and the remaining
nonadditive dimension. Family labels or a high pair-coverage score cannot bypass
this check. Unsupported designs exit 2 without any adapter calls.

`--strategy mechanism_weighted` constrains the old weighting heuristic to those
same four families. It is a post-hoc software repair with no new comparative
validation. It is not the default. The complete negative comparison and its
original software are retained in `research/workflow-qualification-nov01` and
`research/workflow-design-nov02/previous-software`.

Default execution runs four probe pairs per policy, with its own native call:
nine logical queries per policy, plus four shared input-control queries. Five
additional research challenges require `--include-heldout`; two fractional
challenges require `--include-fractional`. With both flags the per-policy union
has 23 logical queries. None of these options executes an adapter without `--run`
or `--qualify`. The earlier 47-query reservation refers to the original four-suite
comparison, whose raw accounting is unchanged.
Use `--native-policy barcode_model` (or your adapter's fixed policy name) to apply
all native controls to a complete calibrated workflow; the compatibility default
is `identity`. A raw uncalibrated decoder need not succeed before a complete
calibrated workflow can qualify. Each scored policy also receives its own native
catalog gate (50 nonblank calls and 99% output support). An unqualified policy has
unsupported held-out truth and its comparative calls are not spent.
Use a new empty output directory for each attempt; failed/partial evidence stays.
The default suite scheduling deadline is three hours, checked between calls;
adapters must bound in-flight execution, including native controls. Each strategy receives four
challenge pairs (eight logical full-workflow calls). Exact cache hits and actual
calibration/decoder executions have separate accounting; equal invocation limits
do not imply identical wall-clock seconds.

Execution produces `summary.txt`, `summary.json` and `recovery.csv`, as well as
the complete control/native/trial JSON. They report native policy, recovery and
composition changes, candidate/reference support, execution failures, unavailable
references, native and fractional gates, and recorded executions/cache reuse.
Results use `observed_failure`, `all_specified_passed` or `unsupported`. Passing
means preserving the declared native catalog on the specified digital challenges;
it does not qualify every possible error or predict natural-error risk.
Per-policy union runtimes are measured over shared queries; separate measured
strategy runtimes are not invented. Failed controls still produce readable
reports and exit status 2. A run with no native-qualified policy also exits 2.
The plan-only command performs zero adapter executions and saves the plan and
readable design audit.

Input NPZ arrays, without pickle:

- `images`, `support`: finite F,Z,Y,X arrays; images>=0, support in[0,1].
- `codebook`: nonnegative F by K array; each column has positive total.
- `rounds`, `colors`: integer vectors F. Two/four-color balanced designs are
  supported, with some additional three-color designs; unavailable designs fail.
- `class_ids`: integer vector K, literal gene/class IDs. Synonymous codewords may
  share an ID. `blank_classes`: integer vector of named blank/control class IDs.
- Optional scalar strings `name`, `preparation`, `metadata_json`. Metadata defines
  ZYX `context`, `calibration_region`, `scoring_region` and `coverage_region` as
  `[start,stop]` pairs. Report and qualify preprocessing before comparisons.

The optional actual BarDensr adapter has correlation and100-iteration regression
backends. Install its pinned upstream dependency in a compatible environment;
planning and the generic API require only the package's normal NumPy/SciPy stack.
Pass `--config config.json` containing `{"decoder":"regression"}` to select it.
The shipped policies are `round`, `frame_template`, `barcode_model`, `group_search`.
They are existing calibration alternatives, not a claim of a superior estimator.

A custom adapter is `--adapter module:factory`. The factory accepts a JSON dict;
its `.run(case, policy)` returns `WorkflowResult(Catalog(points, classes), maps,
coverage, diagnostics)`. Points are array-index ZYX, maps are F by3 corrections
with first-readout gauge fixed BEFORE resampling/decoding. Implementations must
honor the context metadata: the harness adds a16-pixel known-zero XY margin and
updates contexts. No injected maps or reference catalogs are passed. Private
historical-reference provenance paths are removed from the adapter-facing case.
This API contract is not a security sandbox for arbitrary plugin code.

Known-inverse references invoke the SAME complete workflow on the restored input.
Integer images are restored exactly, including their margins. Fractional images
use the actual composed linear resampling, whose loss relative to native is
reported separately; interpolation damage is never called exact inversion.
Empty/insufficient references are unsupported, failed supported predictions are
failures, and missing blanks are unavailable rather than zero error. Object
adapters must explicitly declare candidate support coverage and retain bounding-
box voxel coverage too. Neither conditional support nor blank proxies establishes
whole-field molecular accuracy.

The reference must supply at least 10 nonblank calls before a candidate execution
or support failure can be judged. An unresolved spatial-first matching component
produces an unavailable sensitivity value, independently of label-aware matching.
Integer oracle/native disagreement is an unsupported control outcome. Fractional
oracle-versus-native loss is reported, without subtracting interpolation loss
from the challenged catalog's error. Native count/support gates alone are not a
universal blank-proxy or biological-truth qualification; declare any additional
application-specific quality requirement before comparisons.
