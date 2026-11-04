# Round/color pooling repair — post-hoc development specification

Freeze before observing this repair's results. Run only if the field-one color
diagnostic finds a lower oracle objective and oracle-start identity F1 >= 0.90
for both kernels; otherwise diagnose the objective instead of this initializer.

The image template fitted each of 40 sparse readouts separately, while the
injected errors have shared round or color structure. Pooling observed photons
by these acquisition axes can stabilize initialization without using imposed
offsets, gene identities, barcode assignments, reference calls or outcomes.
This is a conventional structured registration alternative, not a new temporal
measurement or a unique-method superiority claim.

Fit an additive map m[r,c] = a[r] + b[c], with both components centered. Initialize
a from bounded round-sum correlation. Alternate ten round/color updates.
For each update, remove the other component, sum the corresponding axis,
normalize each aggregate by interior L2 energy, build a warped mean template,
and fit bounded correlations on the fixed observed [24,104) interior. Check
support at each step. Keep the existing +/-10 correlation bound. Apply integer
rounding only at decoding, as in the primary study. No outcome-based tuning.

Record two policies separately: this image-only additive fit, and the author's
unchanged 50-step barcode optimizer initialized by that fit. Do not choose between
them using test outcomes or silently select the author loss when it harms outputs.
Run every zero/common/color case on both exposed fields and both frozen decoders:
24 records, plus the original alternatives retained for comparison. Use identical
96-pixel decoding and 80-pixel scoring windows, cuts, support and matching radii.
Reuse calls only after byte-identical full decoder-input verification. Preserve
blank/yield/composition/support metrics and all execution or quality failures.

Development goal: every image-only nonzero case reaches >=0.95 native identity F1
and composition TV <=0.05; every zero case >=0.95 F1. This is not RNA accuracy or
independent transfer. If it fails, retain the result and diagnose before changing
the design. The author-refinement arm is a separate sensitivity, not part of this
pass rule. Synthetic balanced-barcode controls verify translation sign, support,
round/color decomposition, and identity recovery before public repair scoring.

The model cannot represent arbitrary readout-specific/nonrigid residuals.
Two BARseq3 fields remain one exposed preparation. All original primary/sensitivity
results remain; no manuscript superiority or deployment claim follows from this
development repair. A fresh preparation and frozen full workflow are still needed.

Implementation correction before any public repair scoring: the original
unbounded phase initialization violated observed-support safeguards on a
synthetic control. Its exact implementation and failure log are preserved as
structured_repair_v1.py and structured-controls.log. The corrected initializer
uses the same declared +/-10 bounded correlation as every subsequent update.
No public repair outcome informed this change; the full synthetic gate is rerun.
