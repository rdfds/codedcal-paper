# Supplementary Note 2A. Covariance verification and ratio support

## Indices and the count-only counterexample

Panels $(0,1,2,3)$ in Note 2 mean $(a,b,c,d)$ in Note 1. The index $\ell$
in a separable coefficient tensor labels tensor terms, not individual events.
Code-name permutation symmetry acts on the $K$ code labels; permutation of
the three panel arguments is a separate action on three positions.

For the count-only toy, let $b_c$ be the number of unit sources assigned to
code bin $c$, so $N=\sum_{c=1}^K b_c$ is already observed. This is distinct
from the source-panel counts $n_{st}$. At $N=3,K=4$, enumeration of all
$4^3$ assignments gives

\[
\Pr(W=-1)=6/16,\quad\Pr(W=3)=9/16,\quad\Pr(W=27)=1/16,
\qquad E[W]=3,\quad E[W^2]=51.
\]

To make the variance comparison precise, write $C_g=W_g\Delta_g$, with
$(W_g,\Delta_g)$ independent and identically distributed across groups,
$\Delta_g$ independent of $W_g$, and $0<\operatorname{Var}(\Delta_g)<\infty$.
The mean of the individual ratios $C_g/W_g=\Delta_g$ is well defined.
The ratio of pooled sums has first-order asymptotic variance

\[
\frac{E[W^2]}{E[W]^2}\frac{\operatorname{Var}(\Delta)}{G}
=\frac{17}{3}\frac{\operatorname{Var}(\Delta)}{G}.
\]

This is a delta-method comparison, not an assertion of finite-sample ratio
variance. A finite pooled denominator can equal zero (for example, three
values $-1$ and one value $3$); its undefined outcomes must be retained.
The toy also has $E[(W-3)^4]=20832$. The exact law and fourth moment are
checked by rational-arithmetic tests in test_temporal_reporting_models.py.

## Exact specification of the Monte Carlo check

For each exogenous world, the fixed sixteen-group prefix produces an
eight-component vector: $(C_A,R_A,C_B,R_B)$ at each of the two frequencies,
after the implemented placement and angle averages. For each of $M=512$
independent code assignments, let $X_m$ be that group-mean vector and let
$\mu$ be its exact conditional mean. Set $Y_m=X_m-\mu$. The recorded matrix is

\[
\widehat\Sigma_{\rm MC}=M^{-1}\sum_mY_mY_m^{\mathsf T},\qquad
e_F=\frac{\|\widehat\Sigma_{\rm MC}-\Sigma_{\rm exact}\|_F}
{\max(\|\Sigma_{\rm exact}\|_F,10^{-300})}.
\]

It is a known-mean-centered second moment, not a sample covariance centered
at the Monte Carlo mean. For each diagonal component, with
$\overline{Y_j^2}=M^{-1}\sum_mY_{mj}^2$, the estimated Monte Carlo error is

\[
\widehat{\operatorname{SE}}_j=
\sqrt{\frac{\sum_m(Y_{mj}^2-\overline{Y_j^2})^2}{M(M-1)}}.
\]

The standardized difference uses this error with a numerical floor of
$10^{-300}$. The 23 of 864 differences beyond three such errors and maximum
5.38 are retained. Cubic-statistic second moments involve sixth-order source
terms, and these error estimates need not be Gaussian or stable. The toy
above illustrates possible skewness; it does not diagnose the cause of those
23 discrepancies or convert them into passing checks.

Independent exact-enumeration tests compare complex means at relative and
absolute tolerances $2\times10^{-12}$, covariance at relative tolerance
$3\times10^{-12}$ and absolute tolerance $3\times10^{-11}$, and shipped
four-code coefficients at absolute tolerance $3\times10^{-13}$. These are
test acceptance tolerances, distinct from the Monte Carlo error estimates.

## Why an influence-variance gain does not guarantee a usable ratio

For an arithmetic ratio, a necessary condition for a bounded Fieller set is

\[
\overline R^{\,2}>t_{G-1,\,0.975}^{\,2}\,s_R^2/G.
\]

Positive support and the original numerator/reference rules impose additional
conditions. The geometric construction applies the corresponding marginal
condition at its allocated critical value in both windows. The influence
floor $A=\operatorname{Var}(m_C-a m_R)/\mu_R^2$ can cancel common numerator and
denominator variation while the reference floor
$A_R=\operatorname{Var}(m_R)$ remains nonzero. At $B=KG$, its contribution to
the variance of the reference mean is $A_R/G=KA_R/B$. Fewer groups can thus
harm denominator support even when the ratio-influence coding variance falls.
This is an analytical mechanism, not a fitted causal decomposition of the
observed usability changes.

The fixed-state budget experiment samples the same sixteen reconstructed
groups per world for every alphabet size. This deliberate prefix defines a
finite population; it can miss rare groups in the complete acquisition,
especially the difficult EDF10 replicate 2. Its results do not establish
full-acquisition efficiency. The earlier first-pass population scores mean
the per-world `designs` records in results-attempt-01/*.json; the final
finite-state results in finite-budget-attempt-01 and report-final replace
them for reported performance, while retaining those exploratory records.

## Numerical replay across environments

The historical temporal-variant verifier first required relative and absolute
tolerances $10^{-12}$. One already-unusable EDF10 replicate 2 Gaussian
group-common case (replicate 22, $q=0.15$, observed coded assignment 0) differed
by $1.4407731641\times10^{-12}$ relatively between raw and normalized float64
summations; its reference condition measure was 7619.3549. The amended point
check uses the pre-existing arithmetic tolerance (relative $10^{-10}$,
absolute $10^{-12}$). Both the failed check and amendment remain archived;
no uncertainty or outcome gate was changed.

The retained covariance interval recount passes in Python 3.10.16, NumPy
1.23.5 and SciPy 1.10.1. A recount in Python 3.13.2, NumPy 2.5.3 and SciPy
1.18.1 failed an endpoint comparison at relative tolerance $10^{-8}$ and
absolute tolerance $10^{-9}$: maximum absolute difference
$2.79230191\times10^{-7}$, relative difference $2.00041263\times10^{-8}$,
near endpoint 13.95863. Its trace did not record the particular world. A
specific dependency cause has not been established. The failure remains
preserved; the tolerance was not broadened to hide it. These recorded
environment limits are distinct from algebraic validity or physical validation.
