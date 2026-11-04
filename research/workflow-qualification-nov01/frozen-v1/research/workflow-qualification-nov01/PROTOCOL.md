# Workflow qualification by controlled readout challenges

1 November 2026. Development specification, baseline commit
3789b2973740f6660846b951228095e6f9966496. No new comparative outcomes have been
examined when this specification is written. All previously examined BARseq3,
STARfinder and 4T1-tumor preparations remain development evidence. Their old
failures and the registration-certificate experiment are immutable.

## Question and target

At a fixed budget, can a readout-aware computational challenge suite qualify
workflows for preservation of their molecular catalogs on held-out challenges
more usefully than ordinary QC and strong conventional stress testing?
The target is computational consistency under known interventions, not RNA truth,
physical motion, a calibrated natural-error probability or universal safety.
The previous observational risk forecasts are not reused as validated predictors.

A workflow adapter accepts only images, support, codebook/readout metadata and
fixed parameters. It returns a catalog, actual post-correction diagnostics and
execution counts/times. It receives no injected maps or reference catalog.
Full workflows have their own native reference; a separate decoder-only native
reference is retained for qualification and historical comparisons. Empty or
unqualified native catalogs cannot earn a trivial consistency pass.

## Reusable computation contract

Input arrays: nonnegative finite images F,Z,Y,X; matching support; nonempty
nonnegative codewords F,K; round/color vectors F; literal gene/class and blank
mapping. Decoder references and kernel contexts are fixed by the adapter.
All corrections use one reference-readout gauge BEFORE nonlinear decoding.
A scored workflow's known-inverse reference is computed on the known inverse
resampling of the challenged input, in exactly that gauge and decoder context.
This includes unavoidable fractional interpolation rather than subtracting an
optimistic interpolation loss from the outcome. Also retain oracle-versus-native
loss separately. Integer native/no-effect/inverse arrays must match exactly on
full decoder context; fractional controls use analytic interpolation checks.

Reference counts, image support and all adapter gates are qualified BEFORE a
comparative freeze. Native minimum: 50 nonblank calls; scored reference minimum:
10 nonblank calls; at least 99% support on the designated scoring domain.
A missing/failed input control is an adapter failure, not a successful workflow.
Workflow execution failure, support loss or empty predictions with a supported
reference are failures, not silently omitted records. Unavailable blank/QC
measurements remain unavailable and reduce qualification coverage.

Primary consistency thresholds are identity-aware catalog loss >.05 or gene-
composition TV >.05. Match within 3 pixels, independently by gene/class with
maximum-cardinality matching. Also retain spatial-first matching sensitivity,
1.5/5-pixel sensitivity, yields, blanks, support and native self-consistency.
No cell-type or differential-expression significance result is implied.

## Budget and comparisons

Each probing strategy gets FOUR challenge pairs per workflow: challenged and
known-inverse-reference workflow calls (eight logical executions), plus shared
native/control calls. Logical and actual executions, cache reuse, calibration,
decoding and complete wall time are recorded separately. Reuse is never counted
as a new run or assigned inherited runtime. Ordinary no-probe QC is a cheaper
baseline and is labelled as such, rather than described as budget matched.

Strong probing baselines: conventional round/color factorial probes; mixed
round/color/interaction/per-readout probes; fixed-seed balanced per-readout
stress; and geometric/molecular sensitivity summaries from the same queries.
Actual post-correction fit and blank/yield gates accompany the no-probe baseline.
The readout-aware proposal greedily covers codebook- and signal-weighted pairs
of co-occurring readouts using a deterministic balanced-mask bank. It includes
round and color challenges, then chooses two further masks without evaluating
held-out outcomes. This is a standard coverage heuristic applied to this task,
not a new mathematical optimization principle or a calibrated risk model.

Integer translations are the exact-inverse primary test. Fractional translations
are a separate planned panel with oracle resampling references and qualified
subpixel adapters. No claim about natural residual-error prevalence follows
from either panel. Exact seeds, masks, amplitudes, adapter parameters, held-out
challenge bank and all thresholds will be saved and committed before development
comparisons. No outcome-selected amplitude or threshold is permitted.

Development qualification includes exposed BARseq3 fields zero/one and the
already exposed STARfinder culture/tissue inputs if their original declared
conditional targets and native controls can be reproduced. Different fields
of one preparation do not count as independent preparations. Native stock
correlation and nonnegative-regression backends are kept distinct.

The development benchmark is bounded to <=1,600 logical full workflow executions
and <=3 hours of local worker time. No paid/cloud jobs, purchases or new contacts.
If implementation repairs are needed, retain failed controls/outcomes and mark
revised runs as development. Do not search external outcomes for a winner.

## Confirmation gate

Before confirmation, finish the common API, numerical controls, development
comparison and selection rule. Qualify one new publicly accessible preparation
using metadata, native calls/support and inverse controls ONLY. Then commit the
complete comparative freeze before examining its method/challenge outcomes.
The public candidates are Harvard Dataverse 10.7910/DVN/NE0WOP (4T1 culture;
not the previously exposed Zenodo tumor release) and Zenodo 8253771 human M1C
DART-FISH, provided exposure audit and adapter qualification actually pass.
They are candidates, not ready/independent results. Do not open reserved kidney
fields or D076 inputs. If both candidates fail access/qualification, preserve
that limitation while completing all unaffected suite/development work.

A useful result requires improved detection of held-out workflow failures at
comparable false-alert rate, nontrivial acceptance/coverage and declared compute
relative to the strongest conventional suite. Reject-all has no utility. Ties,
losses and unsupported preparations stay explicit. A negative result does not
justify relaxing thresholds or claiming a stronger Bioinformatics paper.

## Implementation freeze clarifications before any comparative outcome

Native v3 corrects the support denominator for the declared object target
(eligible source candidates) while retaining voxel support separately; see
SUPPORT_MEASURES.md. Native tissue's131 correlation and150 regression calls pass
this candidate-domain99% gate. The culture backends remain excluded (candidate
coverage below99%; regression also44<50). BARseq3 correlation inputs have70/125
calls. Image-regression BARseq3 remains explicitly unrun for compute feasibility.
The planned24 development workflows therefore include16 qualified workflows and
8 source-gate exclusions. All six input/backend gates and exclusions are retained.

Compare equal maximum full-workflow invocations, not claim identical observed
seconds. Record actual wall/calibration/decoder time and execution/cache counts.
Primary analysis uses only native-qualified workflows for held-out accuracy;
native workflow failures and every unsupported/unrun/source-excluded case are
reported separately, never silently converted to adequate held-out examples.
Learned native QC includes a training-selected union of up to two observable
scalar gates with>=80% training coverage and the same>=5 adequate-example rule.
The small number of preparations supports descriptive comparisons, not precise
population confidence or independently replicated per-workflow samples.
