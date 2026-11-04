# Annotated calibration decision comparison — 5 October 2026

**A frozen image-derived rule detects the large channel-offset failures on a
held-out annotated crop, across two actual decoder families. It is not yet
independent-preparation validation or a safe deployment recommendation.**

This study uses the BarDensr authors' public hand-curated BARseq spots, not
our decoder outputs as reference. Pinned author correlation and nonnegative
regression kernels both pass unperturbed development qualification. Calibration
is crossed with no correction, common-round phase correlation, a stronger
common-round iterative template, per-frame iterative template, the author's
barcode-model optimizer, and inverse known added offsets. The strong common
baseline was added after development revealed the phase comparator's weakness,
and before opening validation. Its optimizer and ten-step budget match the
per-frame template. All 312 final arms complete; exact call counts and spatial
matching cardinalities independently verify. Failed unbounded development
alignment is retained separately; it is not erased by the repaired 312 arms.

## Frozen result

The rule recommends per-frame rather than strong common-round correction when
the image-derived RMS disagreement of their fitted maps is >=2 pixels. It uses
no imposed displacement or annotation labels at decision time. Development
selects this cut from a declared finite grid; neither thresholds nor the rule
are refit on validation. The target is an identity-F1 improvement >=0.05.

| Validation endpoint, 12 nonzero conditions | Correlation | Nonnegative regression |
|---|---:|---:|
| True beneficial cases detected | 2/2 | 4/4 |
| True negatives | 8/10 | 8/8 |
| False positives / false negatives | 2 / 0 | 0 / 0 |
| Balanced accuracy | 0.90 | 1.00 |
| Mean strong common-round identity F1 | 0.6903 | 0.6572 |
| Mean frozen recommendation identity F1 | 0.7700 | 0.8147 |
| Mean always per-frame identity F1 | 0.7407 | 0.8223 |
| Mean author barcode-model identity F1 | 0.7199 | 0.8153 |
| Recommendations causing >=0.05 harm vs common | 1/12 | 0/12 |

Rule confusion is unchanged at 1.5, 3 and 5-pixel matching radii and under the
label-constrained matching sensitivity. The extra per-frame benefit depends
on decoder: the two 2-pixel component color cases benefit regression but harm
correlation on validation. One false correlation recommendation loses 0.087 F1.
The image feature alone therefore does not certify a universal decoder-independent
benefit. Relative to the author model optimizer, the rule's mean advantage is
0.0500 for correlation and -0.0006 for regression. It does not uniquely defeat
every strong alternative.

The scientific finding is narrower: comparable conventional optimizers can
handle common-round displacement, while channel-specific displacement can
require a richer correction and changes the recovered identity endpoint. The
observation-based disagreement detects the severe cases here. This is a
consequential calibration/decoder interaction, not evidence that channel
registration or cross-correlation is new mathematics. The independently frozen
BARseq3 partial adapter failed native qualification on its second field;
no perturbed molecular comparison ran. A separately frozen conditioning repair
is in progress on that exposed preparation and cannot yet support this finding.

## Controls, scope and remaining weaknesses

Known count-preserving label exchange changes identity agreement while leaving
spatial agreement and composition intact. A known geometric displacement and
inverse correction exercise recovery. Synthetic image controls verify the
author model's correction sign; the previous unbounded template failure is
preserved. Oracle qualification passes on development. Digital resampling can
denoise images and improve validation correlation beyond the native zero
baseline, so inverse-added-offset results are not physical accuracy ceilings.

The reference contains only author 'good' annotations inside [24,176) of each
200-pixel crop: 100 development and 80 validation. Border and other statuses
remain recorded. Unmatched calls may include unannotated real molecules; F1 is
agreement with this incomplete reference, not true molecular false-discovery
rate. No qualified blank-derived FDR is available for these benchmark books.
Spatial ambiguity is recorded, and tighter/label-constrained sensitivities
remain. A twofold threshold change is a stress check; 1.70 for correlation is
outside its attainable evidence range and cannot prove practical robustness.

The 13 conditions per crop are zero and two seeds at component RMS 1,2,4 pixels,
with common-round and common-round-plus-static-color shifts. Many perturbations
on one image are not independent biological replicates. The two crop names do
not establish independent animal/preparation provenance. This is prospective
within-study outcome freezing on retrospective public images, not prospective
instrument validation. No random readout, temporal-budget or cell-type endpoint
was validated here. The prior ISS failure and two temporal modality classes
are unchanged.

Runtime includes both calibration fits to obtain the decision feature; the
rule is not a cheap preflight or a demonstrated speed improvement. Peak extraction,
sensitivities and network are excluded from recorded core timings. Runs used
concurrent local jobs and different declared TensorFlow thread settings, so
stage times are descriptive and do not establish portable speed ratios.

Reproduce: `../decoder-venv.nosync/bin/python run.py qualify`, then development,
`strong_common.py development`, `analyze.py freeze` with the analysis environment,
and commit the freeze before `run_validation.py`. After both validation runners
complete, `../merlin-venv.nosync/bin/python analyze.py report` independently
checks all counts/cardinalities and writes portable call tables; `interpret.py`
produces the complete interpretation and figure. Run from the repository root
with full script paths. Sources are pinned in PROTOCOL.md and source.json;
no external review has occurred.

Disk-preserving projection retains exact portable calls at the primary, half
and twice thresholds, their hashes, density maxima and original tensor hashes.
`compact_caches.py` verifies each projection before removing only its ignored,
regenerable density tensor. `verify_portable.py` recounts all 312 primary arms
at all three matching radii directly from those retained calls and raw author
annotations. The frozen runners regenerate tensors if tensor-level checks in
`analyze.py report` are needed; raw inputs and scientific results remain intact.
