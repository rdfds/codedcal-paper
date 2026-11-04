# Multistart template diagnosis across exposed sources

The 3D frame-template comparator ties the group repair on STARfinder, whereas
initial 2D alternatives failed. Before diagnosis outcomes, test whether this
comes from dimensional information or initialization. Both sources are exposed
and this is post-hoc development, not independent validation.

Use a 2D unnormalized DFT template from three fixed starts: zero, merged-round
reference, first frame. Ten iterations per start; mean-center integer gauge,
reject updates exceeding 10 pixels, choose normalized template residual on
fixed yx16:112 interior. No coding identities, known offsets or reference calls
enter the fit. Qualify exact synthetic zero/common/color translations first.

STARfinder: fit on selected reference plane in crop 0, transfer to crop 4; also
fit on crop 4 itself, matching the 3D comparator observation privilege. Reuse
all three unchanged object decoders and fixed cuts. BARseq3: fit and score each
of both previously exposed fields, using actual prior kernels/cuts and books.
Cases zero, common and color; all outcomes retained. Oracle input/call reuse is
explicit only after exact supported-kernel-input equality. Independent RNA truth
is unavailable. This diagnostic cannot establish a deployable recommendation.
