# Added 3D conventional comparators — exposed development

Freeze this extension after the 54-arm feasibility probe and before extension
outcomes. The positive color result does not include the source 3D comparator.
No threshold, endpoint, scoring window, core or initial probe is changed.

Add three policies to all three cases and three decoder families (27 arms):
source_round_3d: full ROI 3D unnormalized DFT of each merged round against round
one, subtract the mean global translation gauge, then the previously qualified
four-level/10-iteration Demons port on merged round images. Shared map per round.
frame_template_3d: three fixed starts (zero, merged-round registration, and
first-frame registration), up to 10 iterations per start of full ROI 3D unnormalized DFT of each separate
frame to the mean aligned frame template, mean-center translation gauge. Select by normalized template residual on the
fixed z4:26/yx16:112 interior. Reject template starts/updates with >3 z or
>10 y/x displacement. The original symmetric zero-start synthetic test stalled
in a two-mode mean-template basin; retain its failed control log and qualify
this multistart version before native results.
frame_template_local_3d: the same frame template fit plus per-frame four-level
Demons against the fitted mean template. This intentionally gives conventional
methods all scoring-window observations; neither identities nor imposed maps
are available. Local flows are fitted on the scoring ROI rather than transferred
from the calibration ROI. All 30 planes, unchanged complete book/cuts are used.

These are declared Python ROI ports, not an exact MATLAB author reproduction.
Global/local source MATLAB references are retained. Include exact zero/common
translation synthetic controls, frame geometry controls, finite/sign-checked
local flows, and unchanged decoder controls. Preserve any source-local harm.

Report alongside the 54 arms; do not retroactively put these policies into the
original selection rule. A new rule or independent study would require its own
freeze. F1 is recovery of the native conditional catalog under controlled shifts,
not physical RNA truth. This does not establish a natural-error boundary, safety
outside the balanced integer class, or independent preparation transfer.
