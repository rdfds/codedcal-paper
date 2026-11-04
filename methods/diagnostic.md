# Supplementary operational diagnostic

Let $r_A=\overline C_A/\overline R_A$, $r_B=\overline C_B/\overline R_B$
and $D=r_A/r_B$. A diagnostic is invalid unless both sample references and
both window ratios are positive. With $G$ group contributions, set

\[
I_g=\frac{C_{gA}-r_AR_{gA}}{\overline R_A r_A}
-\frac{C_{gB}-r_BR_{gB}}{\overline R_B r_B},\qquad
s_D=\frac{\operatorname{sd}(I_g)}{\sqrt G},\qquad
z_D=\frac{\log D}{s_D}.
\]

The sample standard deviation uses divisor $G-1$. For positive $s_D$, a finite
$z_D$ is quiet when $|z_D|\leq t_{G-1,0.975}$. For $s_D=0$, the implementation
uses $z_D=0$ only if $r_A=r_B$ and otherwise infinity; the latter is invalid.
The arithmetic and geometric stress summaries use this same diagnostic.
A quiet diagnostic is a necessary-model check, not a physical certificate.
The reported quiet-error counts additionally require absolute ratio error
greater than 0.01 relative to the stipulated physical target.

For original/repaired digital comparisons, all source qualification, native
control, proxy and fractional-reference failures remain unchanged. In situ
sequencing means sequencing barcodes at their tissue locations; the sagittal
atlas mechanism test and lag-window sensitivity analyses are separate archived
studies, not additional validation of the temporal motion target.
