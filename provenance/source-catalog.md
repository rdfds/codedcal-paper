# Source Data and analysis provenance

In-tree evidence references, with original records retained and subsequent
analyses explicitly dated. These are project analyses, not external literature
or external review. The ledger indexes 72 entries, including post-hoc composition,
reference reconciliation and conditional temporal validation/variants. Entry 72
is an unperformed protocol retained for provenance, excluded from the submission's
empirical evidence and completion requirements. The detailed older
Results/Discussion remain in `supporting-results.md`; newer development records
below are not automatically new manuscript evidence or independent confirmation.

1. `oracle-crossval-sep15; hardening-sep15`

2. `authors-baselines-sep15`

3. `field-survey-sep16; replication-cd20-sep15; replication-edf1112-sep15; edf10-replication-sep15`

4. `baseline-shaw-sep15`

5. `count-balanced-calibration-sep14; identifiability-theorem-sep16, Thm 1`

6. `d-limit-sep18`

7. `oracle-crossval-sep15`

8. `cleanroom-sep16`

9. `molecular-concordance-sep15; heldout-nup96-sep30` — pooled molecular fit and retrospective pore holdout, including the failed 1% held-out detection rule, released-reconstruction control, train-only diagnostic, and Shaw qualification.

10. `steelman-sdc-sep16`

11. `codedfish-port-sep16`

12. `prospective-blank-sep16; fallback-fix-sep18`

13. `selfaudit-sep19`

14. `novelty-blindcal-sep18; identifiability-theorem-sep16`

15. `certificate-sep18; reliability-audit-sep25`

16. `code-design-sep16; code-tradeoff-sep23`

17. `identifiability-theorem-sep16; prospective-package-sep24; tool/src/codedcal/retrofit.py`

18. `identifiability-theorem-sep16, Addendum 5`

19. `merlin-h2h-zhuang-sep23`

20. `merlin-h2h-sep19`

21. `magnification-estimator-sep19; freeze 1b10c44, with the evaluation script committed before any held-out result existed`

22. `magnification-tms-sep20; frozen before any of those fields was processed`

23. `barseq-crosstech-sep21`

24. `single-round-qc-sep21; freeze eaf434d`

25. `zhuang-release-sep21`

26. `downstream-consequence-sep21; freeze b7075df`

27. `zhuang-release-sep21`

28. `zhuang-release-mouse2-sep21`

29. `zhuang-sagittal-sep24`

30. `zhuang-mechanism-sep25; reliability-audit-sep25`

31. `magnification-estimator-sep19`

32. `allen-audit-sep16`

33. `atlas-scale-sep16; atlas-crossanimal-sep17`

34. `robust-tracker-sep18; wide-window-sep19`

35. `atlas-crossanimal-sep17`

36. `robust-tracker-sep18; wide-window-sep19`

37. `zhuang-wmb-sep18; freeze d736c3d`

38. `released-audit-sep19`

39. `chromatic-artifact-sep19`

40. `mosaic-artifact-sep20; mosaic-artifact-v2-sep20`

41. `archive-census-sep16`

42. `comparator-sep18`

43. `count-balanced-calibration-sep14; stationarity-diagnostic-sep15; phase-preserving-oct24`
    Phase-preserving extension: `tool/src/codedcal/phase_attenuation.py` and
    `tool/tests/test_phase_attenuation.py`; conditional identity and deterministic
    controls, not a new empirical study. Historical results remain unchanged.
    Current reporting guard: `tool/src/codedcal/reporting.py` and
    `tool/tests/test_reporting.py`; 0.5.1 rejects unsupported numerical scales
    around unchanged 0.5.0 numerical modules, with the original source/wheel retained.

44. `oracle-crossval-sep15`

45. `molecular-concordance-sep15`

46. `codedfish-port-sep16; atlas-scale-sep16; atlas-crossanimal-sep17; zhuang-wmb-sep18`

47. `magnification-estimator-sep19`

48. `iss-pilot-sep16`

49. `full-evidence-audit-sep26`

50. `tool/src/codedcal`

51. `shaw-sensitivity-sep26`: explicit default/window lag comparison;
    all four Reinhardt and eight SDC raw/corrected inputs complete; 24 fits,
    verified cloud provenance and exact source-verified Reinhardt repeats.

52. `native-comparison-sep26`: completed native task comparison on all 20
    historical first-round fields (15 accepted, five native fallbacks), with
    held-out and adaptive error-proxy thresholds, achieved target attainment,
    runtime, 55 source-verified barcode tables and full run provenance.

53. `early-comparator-sep27`: same-first-round-input comparison with pinned
    authors' Chromagnon numerical core and actual single-plane application;
    20 fields, raw/spot representations, positive controls, source/raw hashes,
    all numerical return flags and complete-field agreement summaries.

54. `early-transfer-sep27`: unchanged first-round estimator on all 45 existing
    TMS fields, 43 accepted across six experiments; exact offline fit replay,
    per-image input lineage and grouped Allen-derived agreement envelopes.

55. `reference-envelope-sep27`: post-hoc all-round grouped agreement audit,
    preserving original and matched within-round configurations, one unavailable
    revised historical field and explicit limits on coverage interpretation.

56. `augmented-design-sep27`: constructive numerical estimand checks, explicit
    failure of one ideal blur observation under independent plane blur,
    sufficient relative-blur/common-blur models and a rank-deficient
    constant-weight counterexample; design-time algebra, not an acquisition.

57. `bead-reference-sep27`: public four-colour physical bead colocalization
    benchmark, all three overlapping views, spatial calibration/validation
    split, retained failed unchanged transfer and post-hoc Gaussian sensitivity,
    full affine reference, all residuals and exact offline replay.

58. `early-task-comparison-sep27`: frozen downstream comparison of complete
    authors' Chromagnon geometry and first-round spot policy on 20 reused
    Allen fields; single-plane exact decoder, five late-port fallbacks,
    section-held-out radii, actual blank proxies, all count curves, runtime
    stages, raw-range lineage, historical checks and offline count replay.

59. `bead-confirmation-sep28`, `bead-repair-oct04`: independent SSBD 208 attempt with protocol and
    execution frozen before pixels; all twelve full TIFF streams and retained
    projections verified, all 48 channel candidate counts, zero raw-MAD
    diagnosis, explicit post-exposure resource amendment and all 72 unmeasured
    optical endpoints in the unchanged attempt. Detector census complete;
    original optical confirmation incomplete. The exposed Oct04 repair retains
    all 12 stacks, 72 primary comparisons and 216 radius-specific endpoint rows,
    exact detector/localizer replay, hashes, selection losses and independent
    map/residual statistics. Six sparse stacks abstain; denser scale improves
    33/36 pairs over translation and 1/36 over affine. All-gates/all-pairs
    criterion fails; optical characterization is complete as development.

60. `overlap-reference-sep28`, `overlap-resolution-sep28`, `overlap-panel-sep28`:
    public repeated-view physical-observation reference, zero-coordinate metadata
    audit, sparse-support diagnosis, frozen 20-pair extension with all failures,
    fixed early code/full Chromagnon/final native maps, accepted-only and fallback
    accounting, exact range hashes, 60 training-chunk replays and 300 score checks.
    Seventeen references qualify; the new accepted-only superiority criterion fails.

61. `joint-error-model-sep29`: task-specific spatial/temporal error synthesis,
    nine exact rational covariance-identification checks, six positive-semidefinite
    witnesses and a hash-recorded retained-schema/measurement audit. Distinguishes
    relative temporal change from final molecular-separation uncertainty; no new
    empirical accuracy result or joint superiority claim.

62. `landmark-decoding-sep29`: frozen matched-plane joint-slope versus direct
    reciprocal-point/RANSAC landmark comparison on all 20 reused Allen fields;
    seven-plane secondary arm, historical comparators, five policy fallbacks,
    full per-gene/radius curves, section-held-out blank proxies and runtime.
    All 100 historical curves reproduce and 300 new curves replay independently.
    Material benefit fails; source hashes, coordinate-order limitation and all
    field outcomes remain explicit. Development, not native superiority.


63. `native-composition-oct05`: new post hoc gene-composition endpoints on all
    twenty completed native fields and five actual native fallbacks. All 55
    retained export hashes, full-source projection records, run/commit markers,
    prepared codebooks and original readout metadata qualify; 180 historical
    policy totals reproduce. Retains all 500 genes, radial counts, frozen
    1/5/10% cuts, common-cut sensitivity, actual proxies and 360 endpoint rows.
    Independent arithmetic verifies every row; five fresh CSV recounts cover
    4,416,133 rows. A separately frozen secondary amendment uses existing
    donor-selected spot-decoder curves for Chromagnon, joint calibration and
    direct landmarks; all 100 selected endpoints and their word/input hashes
    qualify. Native and spot-decoder scopes remain separate. Same-observation
    composition agreement, not RNA identity,
    equivalence, independent acquisition validation or temporal benefit.

## Subsequent internal development records (outside numbered manuscript claims)

- `root-causes-oct01`: all-field geometry/graph diagnostics, abundance-path
  spectral failures, exact rank-loss replays, and an adverse-field-preserving
  diagnostic figure. Adaptive exposed-data development, not confirmation.
- `joint-calibration-oct01`: donor-excluded joint single/mixture likelihoods,
  retained failed strong-prior sensitivity, shared mode-ambiguity guard,
  connected observational support correction and equal decoder consequence
  tests. Later native donor maps and one-animal scope remain explicit.
- `source-presence-oct01`: explicit inactive-source Gamma-integrated reference
  repair on all 24 retained isolated/absent known-truth controls. Simulation-only,
  one source; no physical identity or multi-source reference qualification.

These records do not add manuscript superiority, molecular or prospective
validation claims. Original numbered Source Data and six main figures remain.

Internal development follow-up (not a numbered primary Source Data item):
`outputs/code-aware-decoding-oct01/RESULTS.md` completes identical code-aware
and conservative incremental decoding across all 20 exposed fields and five
calibration arms. Both 10% practical gates fail. Replacement raises counts
in both arms but fails field error-proxy ceilings; conservative recovery meets
all 20 ceilings with only 0.16% median yield advantage over confidence and
0.26% over separately tuned late port. Independent count/packing audits and
outcome-informed finite-grid diagnostics are retained. Single-plane crop,
blank proxy, adaptive same-animal development; no native, physical-truth or
independent-confirmation upgrade follows. Six main figures remain unchanged.


Subsequent shared-decoder repair (internal development, not a numbered primary
Source Data item): `sparse-photometry-oct01/RESULTS.md` reports all 20 donor-only
photon curves,81.2% median gene-call improvement versus candidate coordinate-only
and 82.1% narrow improvement. All new-arm field proxies≤5%; unique candidate
advantages0.07–0.40% fail every material rule. `dropout-decoding-oct01/RESULTS.md`
and `readout-structure-oct01/RESULTS.md` retain unsafe triplet recovery and modest
warp negatives. All source hashes,25 bands per field,40 donor exclusions,3600
count curves, independent seed/affine and shipped-tool replays, storage attempts
and exact-byte codec/filesystem records remain. Shared established photometry,
not a new temporal observable, native win, true FDR or independent confirmation.
Six main figures and numbered manuscript claims remain unchanged.

Subsequent attribution/headroom audit (internal development):
`calibration-headroom-oct01/RESULTS.md` replays all 3,600 photon curves and
2,880 paired settings. Outcome-informed per-target shared selection gives
3.11% median advantage over narrow; no target reaches 10%. This is a finite-
family diagnostic, not new validation or an impossibility bound. The tested
connected-mixture registration fit does not consume the count-balanced temporal
estimate. Its small decoding advantage does not directly adjudicate the temporal
observable, and its shared photon repair cannot be attributed to that observable.
No manuscript superiority, molecular-reference or confirmation claim is upgraded.

Subsequent downstream calibration-impact evaluation (bounded public follow-up):
`calibration-impact-oct04/RESULTS.md` validates label-independent identity scoring
with count-preserving corruptions and geometric/oracle controls; repairs the
two-color adapter, high-ADU correlation, cosine conditioning and decoder-state
confound; and completes 432 final simulated scenes (6,480 arms) plus four frozen
4T1 fields (480 arms). All call/source hashes and metrics replay. Of 32 planned
strong-decoder external tasks, 28 are supported; low-support native field 1 and
all invalid normalization updates remain. Code-minus-iterative registration
median F1 differences on supported nonzero perturbations are +1.24 pp native and
+0.61 pp Serval, with mixed signs. Large call-reference gains chiefly restore
spatial correspondence and do not establish corrected true tissue identities.
Rich/dose simulation forecasts are identical and score 75% on supported external
tasks; the frozen field-zero lookup scores 100%. One acquisition, one decoding
family, known digital errors and call-derived references; no new independent
RNA truth, biological endpoint or temporal-budget decision benefit. Primary
frozen gains and secondary warm relearning remain separate. Four field-1 members
were retrieved before the amended freeze, with no external scientific outcomes
examined; this exposure and all historical failed models are retained. Six main
figures and numbered manuscript claims remain unchanged.


Subsequent stopped public transfer and repair (internal adaptive development;
not a numbered primary manuscript claim): `calibration-decision-oct04/RESULTS.md`
retains an execution-complete field with invalid endpoints and the interrupted
second field. `calibration-repair-oct04/RESULTS.md` repairs dense label matching,
connected/anchored gains, background prediction and singleton readout folds.
All original 285 project / 60 upstream frozen sources verify. Thirty-two focused
tests pass and all 16 new exposed arms' counts and two reference metrics replay.
Expanded field 200 scoring (twelve 128-square patches) gives 181 reference genes
and two blanks, 4.349% adjusted proxy; field 100 training abstains. The observable
rule is -0.2191 mean F1 versus always-iterative across three imposed nonzero
conditions on this one exposed field. Its code-reference sensitivity fails
proxy qualification. The public JSIT kernel is exactly reproduced under Octave;
a documented standard-FISTA correction passes synthetic controls, but neither
rectified nor published signed public-input preparation yields a qualifying
second-family operating point. The complete native JSIT production workflow
and its published empirical results have not been reproduced or refuted.
All finite negative probes, source hashes, execution failures and extra exposed
three-plane range verifications remain. Four assigned kidney fields remain
unopened. No independent RNA truth/preparation confirmation, reliable deployment
recommendation, temporal-budget validation or manuscript benefit upgrade results.

## 4 October public failure diagnosis follow-up

`outputs/predictor-diagnosis-oct04/RESULTS.md`, `summary.json` and the two frozen
coverage-check records retain all intermediate failures. Corrected readout
responses and iterative fallback remove the exposed0.2191 mean F1 regression
by retaining the comparator; no speedup or unique calibration advantage.
Matched nonnegative SGL plus the existing fixed photometric QC passes the pooled
50-gene/10%-blank development gate:123 genes/1 blank,3.20% adjusted proxy. The
histogram-only check fails:152/21,54.38%. Own-nearest agreement is0.6164 F1 and
123 versus181 gene calls, so superior accuracy is not claimed. Raw public inputs,
source densities/candidates and fixed acceptance rules retain hashes/audits;
37 tests pass and the original285 project/60 upstream freeze verifies. Four
kidney fields remain unopened. No independent preparation confirmation,
physical randomized-readout validation, new modality class or NM upgrade.

## 4 October mixture calling and bounded observable pilot

`outputs/mixture-workflow-oct04/RESULTS.md`, `PROTOCOL.md` and `audit.json`
record a candidate calling-stage repair, not a new manuscript endpoint.
Independent readout fits exposed amplitude shrinkage; fitting-only NNLS refit
repairs36 bright known-source controls. All weak and null results remain.
Public in-sample mixture checking fails451/48 (41.89% proxy); protected isolated
checking gives86/1 (4.58%) and0.4851 own-nearest consistency F1, below the earlier
123/1 and0.6164. No mixture contribution is validated publicly. The debiased
public check is training-only, not full scoring coverage. Two2048-cell-bounded
pilot crops supply observed signal/code/crowding summaries; no target-trained
recommendation or external-preparation confirmation exists. Raw-cache/density/
call and execution-version hashes, all failures,36 tests and285/60 frozen-file
verification remain. Four kidney fields stay unopened. Temporal retrofit,
S3/ISS, two-modality and density-dependent qualifications are unchanged.

## 4 October further public prediction repair

`outputs/public-prediction-repair-oct04/RESULTS.md`, `PROTOCOL.md`,
`fold-design.json`, `full-controls.json` and `audit.json` retain the subsequent
development failure. The later mixture check's multi-bit erasures are replaced
with singleton folds; all individual words remain distinguishable. Residual
candidate expansion repairs the missed weaker word in the two specified
real-codebook mixture controls; all eight bright conditions have identity F1=1.
Weak/null outcomes and the earlier unequal-source failures remain.

Complete on/off evidence, exact alternative mixtures and converged fitting-only
amplitudes are required. The scoring half of the exposed liver training crop
has zero fully supported calls. Its image-SSE ratios to zero prediction are
1.3533 (initial singleton),1.6065 (expanded shortlist),1.6685 (wider shape) and
1.9030 (wider shape/trained response/color nuisance). Better calibration shape
fit does not transfer. Fixed inherited geometry/normalizers and all-readout
full proposals make these conditional diagnostics, not end-to-end independent
validation or RNA truth. The newer model was not advanced to the twelve scoring
crops after training failed.

42 regression tests pass;136 prediction hashes, complete support/counts,
source versions and the285/60-file original freeze verify. The final nuisance
model has retained ordered components and a pixel archive that exactly matches
the original22 predictions. Four kidney fields remain unopened. No qualified
workflow recommendation, independent preparation result, physical temporal
measurement, modality expansion or Nature Methods claim upgrade follows.

## 5 October conditional public prediction repair

`outputs/public-prediction-model-oct05/RESULTS.md`, `PROTOCOL.md`, `audit.json`
and the retained fold/penalty/control records govern the latest exposed
development. A separate conventional ridge forecast reduces scoring image SSE
to 0.778719 times zero, with all 22 readouts improved. The old expanded source
feature adds no benefit. Count/composition-preserving label corruption leaves
that image forecast unchanged while identity F1 drops to 0.0811, so the image
endpoint cannot certify molecular output or select a workflow by itself.

GLS nuisance conditioning, measured empirical spatial dispersion, selected-
dictionary profile information and response-balanced candidate pursuit repair
specified known-source failures. Calibration-only sparsity minimization is
constrained by recovery: penalties 12, 24 and 48 fail; six passes all eight
bright conditions and preserves weak/null abstention. The final conditional
source model gives five scoring genes/no blanks under all 22 vector checks,
but image SSE/zero remains 1.806830. The 130 proposals partition into 72 missing
response support, 37 incomplete sources, 16 vector failures and five passes.
Calibration-scaled source increments give 0.800460 versus 0.800334 without them;
that diagnostic supplies no material source-specific image gain. Empirical
correlations can contain signal as well as noise. Only 15 response channels
have model support from 25 proposed isolated sources; seven neutral priors do
not validate a call. Sources, labels and gains are not physical RNA truth.

Seventy-three regression tests pass; 244 new forecasts, 88 historical exact
reconstructions, source/count evidence, retained public folds, penalty selection,
geometric controls and the original 285 project/60 upstream freeze verify.
Default forecasts remain bitwise unchanged in the specified regression control.
All startup, control, source and partial failures remain. Only the previously
exposed liver field-200 training crop was used, with inherited iterative
geometry/normalization and separate calibration/scoring halves. No new public
acquisition, paid compute or reserved kidney field was accessed; the twelve
earlier scoring crops were not advanced after source qualification failed.
No usable action model, independent preparation result, source superiority or
new temporal observable is established. Main manuscript/primary figures and
S3/ISS, modality, density and simulation-only claim limits remain unchanged.

64. **Annotated calibration/decoder decision comparison (5 October 2026).**
`outputs/independent-comparison-oct05/PROTOCOL.md`, settings-freeze.json,
preflight-freeze.json, STRONG_COMMON_AMENDMENT.md, comparison-summary.json,
all-arms.csv, interpretation.json, density-and-call-manifest.json and RESULTS.md.
Source: Chen et al., BarDensr (2021), DOI 10.1371/journal.pcbi.1008256;
author repository commit dfdafe9055be322b286369d21e3349204c1f6fc3 and two
checksum-verified annotated HDF5 crops. Actual correlation and nonnegative
regression kernels cross six calibration policies; all 312 final arms complete.
Independent global assignment verifies every primary spatial cardinality and
call count; portable exact call CSVs and density-cache hashes are retained.
Failed unbounded development alignment is preserved. Before validation, the
strong common-round iterative baseline and 2-pixel observed map-disagreement
rule were frozen. On twelve nonzero validation conditions, correlation gives
TP/TN/FP/FN=2/8/2/0 and regression 4/8/0/0; radius and label-constrained
sensitivities agree. One correlation recommendation harms F1 by 0.087. Mean
selected/common F1 is 0.770/0.690 and 0.815/0.657; the author model gives
0.720/0.815. No unique superiority across families follows. Author good
annotations are incomplete; F1 is reference agreement, not molecular FDR.
Digital offsets, 100/80 interior reference points and unestablished independence
of the two crops limit inference. Full preflight fits both policies and timings
are uncontrolled. Supplementary Figure 4 displays every condition. The separate
BARseq3 resource audit and native qualification are development/unfinished
transfer, not another validated preparation or temporal modality.

The separately specified `INTEGER_CONTROL.md`, integer_control.py,
integer_analysis.py and integer-control/ retain all 120 integer-grid sensitivity
arms with no failures. Native normalization and a common decoder window make
oracle calls exactly equal zero/no-correction calls in both families/crops.
Scoring-crop selected/common mean F1 is 0.7083/0.6363 and 0.8235/0.6137;
author-model values are 0.6516/0.8326. One harmful correlation recommendation
persists. This is post hoc exposed-data sensitivity, not new preparation
confirmation; quantization and the common window are explicitly declared.

65. **External matched-distribution calibration/decoder comparison (6 October 2026).**
`outputs/barseq3-native-repair-oct05/PROTOCOL.md`, RECOVERY_PROTOCOL.md,
REPRODUCTION.md, READOUT_STRUCTURE_ANALYSIS.md, PRIOR_ART_BOUNDARIES.md,
full-input-verification.json, environment.json, decoder-environment.json,
field-0/, field-1/, recovery/summary.json, recovery/all-arms.json,
recovery/RESULTS.md, recovery/readout-stratification.json,
recovery/runtime-accounting.json, verification.json and
independent-correlation-verification.json. Qi et al. BARseq3 preprint,
DOI 10.64898/2026.05.13.724900; Zenodo 19024678; author source pinned at
32df2c4793b8ab0534b73cc56a506ae16d4048a8, BarDensr at
dfdafe9055be322b286369d21e3349204c1f6fc3. Whole-file published MD5 verifies.
Source indices are retained because their FOV-name mapping is unverified.
Native full-field conditioning repairs the earlier partial adapter without
changing decoder cuts/gates. One rather than five unreleased fields supplies
the fixed normalization profile. Native four gates pass with zero observed
blanks and full support, an operational qualification rather than true FDR.

All sixty pre-specified arms complete. Five policies, two actual decoder
families, zero/common-round/color cases and two fields remain. The complete
nonzero offset multisets are identical (twenty +(4,2), twenty -(4,2)); allocation
alone differs. Strong common correction restores both fields' common case
exactly. Color-case correlation F1 is 0.000/0.030 and regression 0.006/0.008;
author model gives 0.908/0.628 on field zero, 0.059/0.012 on field one. The
unchanged two-pixel selector misses field one's color case. That field's strong
round-only regression has 437 nonblank calls, zero blanks and F1 0.008 against
66 native algorithm reference calls. Equally optimal global matching gives
0.01193; label-constrained five-pixel matching gives 0.13519 and composition
TV is 0.87643. Spatial ambiguity affects exact identity numerators; the severe
output change remains. All point-proxy, output-count and support
misses remain. Zero selected harms does not establish successful recovery.
Reference F1 is native-output consistency, not independently known RNA truth;
one external preparation, an exposed adapter field and a spatial holdout do
not supply two preparations or prospective intervention validation.

Exact integer inverse controls reproduce native calls; cached solver-input
reuse is explicit. Retained-call metrics at three radii, twenty native array
hashes and source lineage verify; an independent NumPy correlation/peak kernel
reproduces both native reference lists. Main Figure 4b shows every
condition. FIT_SUPPORT_SENSITIVITY.md separately declares a post-hoc interior-fit
check of artificial borders on the same sixty-arm design; its complete results
and verification are retained in fit-support/RESULTS.md and verification.json and do not replace the primary
cohort. All sixty sensitivity arms complete; development color-template F1
improves to 0.588/0.371, but spatial-holdout values remain 0.015/0.012 and
the frozen selector still misses it. Native hash/recount checks and 42 exact
input reuses verify. Global matching verifies cardinality/minimum distance
for all 120 primary/follow-up arms; 28/360 radius-specific comparisons have
equal-cost identity differences, maximum observed F1 change 0.03361.
CONCLUSION.md synthesizes all outcomes and inference limits. The isolated-source Cauchy bound requires at least seven supported
readouts at the actual correlation cut; 1,151/1,770 code classes have at most six
readouts in either separated camera group. PSF tails/overlap/background prevent
an unrestricted pixel-level prediction. COUNTEREXAMPLE.md retains a selected
actual-kernel absent-class call with no blank call, not a prevalence estimate.
Sparse blank uncertainty and the failed specificity archive schema audit remain
explicit. No randomized temporal intervention or measured T comparison follows.

66. **Post-hoc registration diagnosis and bounded global repair (6 October 2026).**
`outputs/barseq3-native-repair-oct05/OBJECTIVE_DIAGNOSIS.md`,
objective_diagnosis.py, objective-diagnosis/, STRUCTURED_REPAIR.md,
structured_repair.py, structured-repair/, GLOBAL_REPAIR.md, global_repair.py,
global-repair/, REPAIR_RESULTS.md, REPAIR_SUMMARY.json, verify_repairs.py,
repair-verification.json, interpret_repairs.py and repair-comparison.png/pdf.
Same checksum-verified BARseq3 source and actual kernels as Source Data 65.
All six oracle-start cases (twelve decoder records), twenty-four pooling records
and twenty-four global/refinement records complete. The field-one color inverse
has author loss 0.385586 versus 0.496903 at the centered failed geometry;
oracle-start outputs are exact. A nine-point profile exposes a padded-objective
barrier, supporting a local-basin limitation without proving an overlap mechanism.
The pooled initializer still fails field one. Original unbounded-phase and
round-first global-control failures are preserved with exact source versions.

The corrected search compares zero-round and bounded round-fitted families,
all three balanced color partitions, 81 coarse and nine refinement evaluations
per family/partition (540 entries). Photon-preserving composition and a shared
input denominator remove an intermediate-cropping confound. Selected losses
match actual author TensorFlow and independent direct float64 projections.
All twenty-two corrected synthetic controls recover images exactly. Search inputs
exclude imposed maps, partition labels and reference calls. All twelve global
primary field/case/decoder combinations restore byte-identical native inputs,
exact call lists, F1 1.000 and composition TV 0, retaining original yields,
zero observed blanks and full scoring support. Reuse is explicit, not independent
solver replication. Optional local refinement damages field zero's common case
(correlation/regression F1 0.9787/0.8114); all outcomes remain. Retained-call and
whole-domain assignment verification covers all sixty new arms at three radii,
twenty array hashes and six direct global losses. Tie sensitivities remain.

Main Figure 4c displays every repair condition; the original diagnostic
figure remains in this source package.
This is exposed development on one preparation, integer residuals and a bounded
balanced-group model; the exposed amplitude is on the coarse grid. It does not
establish independent transfer, arbitrary per-readout/nonrigid coverage, physical
RNA accuracy, temporal T attribution or unique early-method superiority.
Search runtimes are retained; optional-refinement timing is unavailable and
no controlled speed ratio is inferred. `outputs/barseq-public-transfer-audit-oct06/`
retains new public-preparation metadata, sample-ID inconsistencies and download
timeouts; no candidate intensity data or independent result has been examined.


67. **Public STARfinder qualification, strong-comparator and source-architecture diagnosis.**
`outputs/starfinder-public-confirmation-oct06/` retains complete member CRC/SHA
provenance, source mapping, every failed original native interface, conditional
supported-object cuts/abstentions and actual pinned correlation/nonnegative
kernels. Culture search recovery ties strong 2D/3D image multistart. Independent
recount verifies all 111 initial/added-comparator/initialization rows. Stronger
barcode starts add twelve BARseq3 rows: field-zero color F1 .971/.979, field-one
.015/0; the older .0294/.0086 author policy at radius three remains stronger in field one. Global
matching changes three of 36 new barcode-model radius scores by at most .00709.

All seventeen tissue source members are retained. A complete pre-pixel frozen
workflow fails native eligibility: 18/47/91 calls and severely unsupported maps;
its 96-arm comparison is unrun. Post-hoc full-field source registration repairs
large round shifts before extracting ROIs, rather than fitting them inside small
crops. Its conditional native catalog has 544/668/748 calls, minimum support
.9681; original whole-interface .99 remains failed. All 96 repaired comparison
arms complete. Search and oracle recover zero/common/color catalogs exactly;
strongest currently tested 2D author color F1 is .926/.920/.913. Outside its
balanced-group class, search gives .467/.532/.527 versus author .984/.954/.946.
A separately qualified, committed full-volume author comparator adds three arms
and recovers .979/.958/.956 F1. The search gain is 2.1–4.4 percentage points,
below the specified five-point distinct-benefit criterion. All 297 global tissue
assignments agree with primary F1. The original 96-arm report remains intact.

The warped-pixel decision rule harms zero/common cases and misses color benefit.
Direct unchanged-input object scoring fixes those harms but still misses benefit;
a first-source-color constraint also fails. The latter uncommitted development
pilot's staging dependency error is explicitly recorded. No threshold was lowered
to force a favorable recommendation. All failures, support, identity/composition,
yields, abstentions, negative-code limitations and runtime are retained.
Frozen source libraries and the group-search core are unchanged. All follow-ups
are exposed development, not independent complete-workflow confirmation, natural
RNA truth, a validated recommendation or physical randomized-readout validation.


68. `output-risk-oct06` — experimental molecular-output risk assessment.

Original finite-profile study completes thirteen exposed cases / twenty decoder
records and fails decision qualification. Its narrow assessment eliminates all
culture objects; bootstrap-weighted repair also fails. A separately frozen
common-canvas/full-domain repair passes all twenty actual inverse controls,
retains 860 candidate catalogues and gives 13/13 correlation decision coverage.
Maximum/weighted rules detect 2/2 versus 1/2 failures with 5/11 versus 2/11 false
alerts. Both failures are in BARseq3; no added value over fit-gap/geometric
checks is established. Incumbent-selection diagnosis identifies an older
weaker model replacing a stronger retained author workflow. A label-preserved
author sensitivity completes all twenty forecasts but also fails, with an
explicit severe-output-loss abstention and an unrecomputed confidence comparator.
All original records and failed gates remain. Independent recount verifies
profiles, bootstrap/tie weights, forty common-domain forecast statistics, native
reference assignments and forty preserved-workflow statistics. Conditional
native-output consistency, not RNA accuracy, prospective transfer, a risk
certificate or a randomized temporal measurement. See RESULTS.md, protocols,
freezes, common-domain/verification.json and preserved-workflow/verification.json.

69. `reference-reconciliation-oct07; reconciliation-sensitivity-oct23` — complete
existing coded-acquisition reference audit and full-cohort sensitivity.

Eight distinct acquisitions / nine raw-corrected views include NPC, SDC,
EDF10 replicates 1–2, Figure4 untreated and RTX replicates 1–3. All four
Figure4 RTX member hashes equal RTX replicate 1, so the alias counts once.
Historical source contrasts, eligibility and event weights replay before
interpretation; original code moments are audited with their original phase
precision. Centered float64 diagnostics retain numerical differences separately.

Exact event/localization counts, source selection, temporal panels and source
mean event frames are retained. Assignment-assisted diagonal Fourier reference
and code use identical original selected localization rows. This removes
specific target/population mismatches but conditions on oracle selection;
code exogeneity and physical source identities are not established. Original
grouped-code populations, every historical fit variant, unsupported ratios,
paired spatial-block intervals and all disagreements remain. A bounded fit
or overlapping interval is not an equivalence certificate.

Six CD20 acquisitions supply 96 primary field-derived control worlds / 384,000
group draws, with two noise models, two seeds and four injected variances.
Original source QC is unchanged. All candidate clusters precede source QC;
inferred emitters may already be merged. Declared response stability failures
are retained, as are the initial 400-draw pilot and pre-control import failure.
Paired simulation influences retain common random draws across injected T.
Sparse clustering passes strict-cutoff partition and retained-world checks.

Source selection narrows several old differences; EDF10 replicate 2 retains
a residual difference at both frequencies; corrected SDC differs at q=0.10. SDC's prior intercept removes modeled physical settling,
so its corrected injected-component target is not total physical code motion.
No pooled universal lag, independent physical truth, quiet-D certificate,
unselected-population validation or randomized MERFISH intervention follows.
See RESULTS.md, ESTIMANDS.md, comparison.json, verification.json and
Supplementary Table 2 (`reference-comparison.md`).

70. `temporal-validation-oct09` — frozen preparation-derived temporal estimator
accuracy/uncertainty study; all 126 worlds, 3,528 simulated field realizations and 7,056 coded
observations retained. All pre-QC Source Data 69 candidate clusters, source
catalogs, burst multiplicities, stationary Gaussian/residual laws, physical vs
observable targets, stress failures, finite/unbounded Fieller sets and paired
field uncertainty remain. Source labels/laws are synthetic assumptions, not
physical validation. Complete summary, group records, independent recount,
weight-concentration audit and post-hoc resolution/boundary sensitivities;
see RESULTS.md and Supplementary Table 3.

71. `temporal-variants-oct17` — separately specified retained-observation
comparison of existing arithmetic and geometric temporal constructions. All
126 worlds, 84,672 variant cases and 2,016 condition summaries retained, with
42,336 original arithmetic case/1,008 summary replays and 84,672 independent
geometric-window formula checks. Conservative joint confidence sets, changed
heterogeneous targets, resolution losses, all kernel/code failures and the
precision-only verification amendment remain. Conditional computation, not
new acquisition or physical validation; RESULTS.md and Supplementary Table 4.
The `geometric-integration-oct18` software audit exactly replays all 84,672
retained cases through codedcal 0.4, with analytic/CLI tests and an isolated
wheel check; it adds no scientific observation or automatic variant selector.

72. `public-confirmation-oct17` — complete bounded D076 comparative protocol,
executor, independent recount and artificial controls, frozen before either
preparation's pixels. Two preparations, six fields, seven conditions, six
policies and two actual decoder families; all 504 slots retained. Source
access and real native qualification are pending. This is an executable plan,
not an independent result, a physical temporal intervention or RNA truth.
