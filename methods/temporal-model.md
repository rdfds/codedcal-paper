# Supplementary Note 1. Calibration observables and their identifiability

This current exposition makes the panel averaging explicit and typesets the
same conditional identities as the retained historical derivation. The
historical proof, all counterexamples and numerical records remain unchanged.
Hadamard extraction, cumulants and Fourier dynamical normalization have
established precedents [17-19,35,36]; the increment is the specified
count-product-matched observable under permanent source-code randomization.

## 1. Permanent-code source-diagonal extraction

Source $s$ has position $\theta_s$ and an independent uniform permanent code
$c_s\in\{1,2,3,4\}$. Let $h_1,h_2,h_3$ be the nonconstant four-point Hadamard
rows, with pointwise product $h_1h_2h_3=1$. For arbitrary source contributions
$a_s,b_s,d_s$, define $A_j=\sum_s h_j(c_s)a_s$ and analogously $B_j,D_j$.

\[
H_3(a,b,d)=\frac16\sum_{(j,k,l)\in\operatorname{Perm}(1,2,3)}A_jB_kD_l.
\]

**Lemma 1.** Conditional on contributions independent of the assigned codes,
$E_cH_3=\sum_s a_sb_sd_s$.

**Proof.** Expand over three source indices. Three distinct indices give a
product of independent zero-mean signs. A term with two equal indices has a
remaining zero-mean independent sign. For three equal indices the product
of the distinct Hadamard rows is one. Each permutation therefore retains
exactly the same source diagonal. This is an expectation identity, not exact
extraction at a realized finite code assignment.

Every event from one persistent source must carry that source's same code.
Randomly splitting events after acquisition does not implement this model.

## 2. Basic three-panel matched-count identity

For angular spatial frequency $q$ and disjoint event panels $a,b$ in window A
and $c$ in window B, write

\[
X_{st\ell}=\theta_s+D_{s,w(t)}+\epsilon_{st\ell},\qquad
Z_{st}(q)=\sum_\ell e^{\mathrm{i}qX_{st\ell}}.
\]

The unit-weight count is $n_{st}$. The following conditions are sufficient:

1. Permanent source codes are independent of positions, counts, displacements and event errors.
2. Conditional on counts, source nuisance and persistent window displacements, errors in disjoint panels factor with a stable source-specific characteristic function $\psi_s(q)$ independent of displacement realizations.
3. Conditional on the complete profile determining count-product weights and kernels, $E[e^{\mathrm{i}q(D_{sA}-D_{sB})}]=a(q)$ is the same real positive attenuation for each contributing source.
4. The population reference is nonzero. Sample ratios additionally need their stated sign and numerical-support checks.

The basic population moments, with expectation over codes and observations,
are

\[
\begin{aligned}
R_{A;c}&=E\operatorname{Re}H_3(Z_a,\overline{Z_b},n_c),\\
C_{A;c}&=\tfrac12E\operatorname{Re}\{H_3(Z_a,n_b,\overline{Z_c})+
H_3(n_a,Z_b,\overline{Z_c})\}.
\end{aligned}
\]

**Theorem 2.** Under these conditions,

\[
R_{A;c}=\sum_s E[n_{sa}n_{sb}n_{sc}|\psi_s(q)|^2],\qquad
C_{A;c}=a(q)R_{A;c}.
\]

**Proof.** Lemma 1 removes cross-source products. Position and within-window
displacement cancel from the reference; conditional panel factorization leaves
$|\psi_s(q)|^2$. In either cross term position cancels, the same error factor
remains, and the conditional displacement expectation contributes $a(q)$.
All terms carry the identical three-panel count product. Summation proves the
identity. For fixed kernels its weight is $E[n_{sa}n_{sb}n_{sc}]$; sources
absent from this product are outside the target's support.

For separately normalized pair moments, put

\[
u_s=\frac{n_{sa}n_{sb}}{\sum_jn_{ja}n_{jb}},\quad
v_s=\frac{n_{sa}n_{sc}}{\sum_jn_{ja}n_{jc}},\quad k_s=|\psi_s(q)|^2.
\]

With positive total pair counts,

\[
\frac{C_2}{R_2}=a(q)\frac{\sum_s v_sk_s}{\sum_s u_sk_s}.
\]

The nuisance factor equals one for every admissible positive kernel vector if
and only if $u=v$. Sufficiency is cancellation. Necessity follows because
$\sum_s(v_s-u_s)k_s=0$ for every positive vector; varying one coordinate
forces its coefficient to vanish. Independent Gaussian errors already allow
these variations at fixed nonzero $q$. This is a statement about this pair-ratio
family, not a universal necessity of random codes or third-order measurements.
Known or homogeneous kernels and other observation models allow other methods.
The calculation is not an empirical failure test of kICS, DDM or SOFI. In
particular, heterogeneous-population kICS already removes shared image/PSF
factors with its normalized spectral autocorrelation and a stated blinking/
transport mixture [34], equations 14–15.

## 3. Actual four-panel implementation and pooling

The software uses $a,b$ in A and $c,d$ in B. For $j\in\{c,d\}$ let
$R_{A;j},C_{A;j}$ be the basic formulas with $c$ replaced by $j$. Define
$R_{B;k},C_{B;k}$ for $k\in\{a,b\}$ by swapping the windows, with $c,d$
the within-window phase panels. The implemented moments are

\[
\begin{aligned}
R_A^{(4)}&=\tfrac12(R_{A;c}+R_{A;d}),&
C_A^{(4)}&=\tfrac12(C_{A;c}+C_{A;d}),\\
R_B^{(4)}&=\tfrac12(R_{B;a}+R_{B;b}),&
C_B^{(4)}&=\tfrac12(C_{B;a}+C_{B;b}).
\end{aligned}
\]

Thus the source weights are $n_{sa}n_{sb}(n_{sc}+n_{sd})/2$ in A and
$n_{sc}n_{sd}(n_{sa}+n_{sb})/2$ in B. Under Theorem 2 each cross moment is
$a(q)$ times its own reference despite different weights between windows.
Hats below denote means of the actual group contributions, with the specified
orientation averaging performed before pooling. The two estimators are

\[
\widehat a_{\rm arith}=
\frac{\widehat C_A^{(4)}+\widehat C_B^{(4)}}
{\widehat R_A^{(4)}+\widehat R_B^{(4)}},\qquad
\widehat a_{\rm geom}=\sqrt{
\frac{\widehat C_A^{(4)}}{\widehat R_A^{(4)}}
\frac{\widehat C_B^{(4)}}{\widehat R_B^{(4)}}}.
\]

Arithmetic is a ratio of pooled moments, not a mean of window ratios.
Geometric support requires both original positive usable window ratios and
its original marginal confidence-set rules. Subsequent notation $C_A,R_A,C_B,R_B$
refers to these four-panel moments unless a basic panel triple is explicit.
Finite-sample ratios are not unbiased. Under an independent zero-mean Gaussian
window model with common and source-specific variances $\gamma^2,\tau^2$,

\[
a(q)=e^{-q^2T},\qquad T=\gamma^2+\tau^2=-\log a(q)/q^2.
\]

Finite variance or independence alone does not justify this exponential.
Without Gaussianity, attenuation at the observed frequencies is the target;
a variance interpretation through curvature at zero needs an appropriate
second moment and justified limiting estimator. Internal groups cannot
independently resample a shift shared by an entire acquisition or identify its
across-acquisition variance from one window pair.

## 3.1 What the split diagnostic does and does not prove

Put $r_A=C_A/R_A$, $r_B=C_B/R_B$ and $D=r_A/r_B$. Theorem 2 gives $D=1$.
Rejection is evidence against the model; non-rejection is not a certificate.
The retained population counterexample has $D=1$ but frequency-dependent
apparent variances of 7.8 and 16.0 $\mathrm{nm}^2$ despite zero physical motion;
a separate stability-violation witness exceeds 40 $\mathrm{nm}^2$ bias.
Frequency and orientation checks are necessary diagnostics, not sufficient
physical-validation conditions.

For a source-uniform change $\psi_{sB}=c(q)\psi_{sA}$, the real contrasts obey

\[
r_A=\operatorname{Re}(c)a,\qquad
r_B=\frac{\operatorname{Re}(c)}{|c|^2}a.
\]

If $c$ is positive real, $\sqrt{r_Ar_B}=a$ although $D$ can reject. For complex
$c$, the geometric target is $a|\operatorname{Re}(c)|/|c|$ and generally
fails; a kernel mean shift supplies a counterexample. Positive ratios are
required for logarithmic estimates.

## 3.2 Temporal support and the qualified fourth-order alternative

The four panels do not define one instantaneous lag. For source $s$ let
$W_s=\prod_{t=1}^4n_{st}$ and $\mu_{st}$ be its mean retained event coordinate.
A directly defined signed target is

\[
B_{\rm obs}=\frac{\sum_sW_s(\mu_{s1}-\mu_{s4})(\mu_{s2}-\mu_{s3})}
{2\sum_sW_s}.
\]

Conditionally zero-mean errors with zero cross-panel covariance cancel from
its expectation even when their variances differ. The latent quantity is
half the weighted expectation of $(D_{s1}-D_{s4})(D_{s2}-D_{s3})$; it is half
a covariance only for centered increments. Actual event timing, selection and
count conditioning are part of the model; panel midpoints generally change it.

The four-code kernel $h_4$ has coefficient $1,-1/3,1/3,-1$ when the label tuple
has respectively one, two, three or four distinct labels. For an index tuple
with $m$ distinct sources the label distribution depends only on $m$, not
multiplicities. In particular, source partitions $(3,1)$ and $(2,2)$ both give
probabilities $1/4,3/4$ for one or two labels, so their kernel mean is zero.
Enumerating the $4^m$ assignments for $m=1,2,3,4$ gives means $1,0,0,0$.
The corresponding $H_4$ therefore extracts the fourth source diagonal in
expectation. Exact rational checks retain all five source-index partitions.
For code-bin counts $N_t$ and coordinate sums $S_t$, define

\[
\begin{aligned}
D_0&=H_4(N_1,N_2,N_3,N_4),\\
N_0&=\tfrac12\{H_4(S_1,S_2,N_3,N_4)+H_4(N_1,N_2,S_3,S_4)\\
&\hspace{20mm}-H_4(S_1,N_2,S_3,N_4)-H_4(N_1,S_2,N_3,S_4)\}.
\end{aligned}
\]

Then $E_cN_0/E_cD_0=B_{\rm obs}$, not an unbiased realized ratio.
For phase sums $Z_t(q)$ set

\[
\begin{aligned}
C_4(q)&=E\operatorname{Re}H_4(Z_1,Z_2,\overline{Z_3},\overline{Z_4}),\\
R_4(q)&=E\operatorname{Re}H_4(Z_1,\overline{Z_2},Z_3,\overline{Z_4}).
\end{aligned}
\]

For $u=X_1-X_4$ and $v=X_2-X_3$ the phase patterns are $u+v$ and $u-v$,
whose squared difference is $4uv$. With finite second moments and nonzero
count mass,

\[
\lim_{q\to0}\frac{-\log[C_4(q)/R_4(q)]}{4q^2}=B_{\rm obs}.
\]

At finite frequency, independent symmetric event errors give the same real
product of four characteristic factors in both patterns. A centered jointly
Gaussian latent four-coordinate law with common conditional $B$ then gives
$C_4/R_4=e^{-4q^2B}$. This does not cover arbitrary asymmetric kernels,
changing means or cross-panel dependence and supplies no efficiency guarantee.

For common centered covariance $K(s,t)$ and panel-time densities $p_t$,
independent of source, let $\overline K_{uv}=\iint K(s,t)p_u(s)p_v(t)\,ds\,dt$.

\[
B=\tfrac12(\overline K_{12}+\overline K_{34}-\overline K_{13}-\overline K_{24}).
\]

Independent constant windows recover $T$. Brownian increments with variance
$2\kappa$ per unit time give $B=\kappa(E t_3-E t_2)$ for ordered nonoverlapping
panels. Deterministic linear drift can give positive $B$ with zero stochastic
variance; oscillatory covariance can give negative $B$. Unrestricted
panel-dependent error means can absorb displacements without changing any
observations, and correlated errors can contribute to the same product.
No universal resolution or molecular-distance error bar follows. All original
unbounded/disconnected outcomes and negative controls remain; this alternative
is not promoted to a physically validated production estimator.

## 3.3 Phase-preserving source-uniform complex extension

Keep code exogeneity and conditional panel factorization, with within-window
stable kernels $\psi_{sB}(q)=c(q)\psi_{sA}(q)$ for the same nonzero complex $c$
across sources and groups. Let $m(q)=E[e^{\mathrm{i}q(D_{sA}-D_{sB})}]$ be common
conditional on count/kernel weights; it need not be real. Require nonzero real
references. Tildes denote complex cross moments with the same four-panel
placements, before taking real parts.

**Corollary 2.1.**

\[
\frac{\widetilde C_A}{R_A}=m\overline c,\qquad
\frac{\widetilde C_B}{R_B}=\frac{\overline m}{\overline c},\qquad
\frac{\widetilde C_A\widetilde C_B}{R_AR_B}=|m|^2.
\]

Thus $\sqrt{|(\widetilde C_A/R_A)(\widetilde C_B/R_B)|}=|m|$.
Source-diagonal extraction is complex multilinear. In A, cross and reference
error factors are $\overline c|\psi_{sA}|^2$ and $|\psi_{sA}|^2$ with the same
weights. In B they are $c|\psi_{sA}|^2$ and $|c|^2|\psi_{sA}|^2$, with reversed
displacement factor $\overline m$. Division and multiplication prove the claim.

The retained 5-nm kernel-shift witness with $T=9$ has real moments at $q=.10$
indistinguishable from a stationary Gaussian model with $T=22.058\ldots$.
Keeping complex moments recovers 9 at both retained frequencies instead of
22.058/22.884. This exact population/enumeration control is not another
acquisition or a sampling-coverage study. For Gaussian increments with mean
$\mu$ and variance $2T$, $m=e^{\mathrm{i}q\mu-q^2T}$, so $-\log|m|/q^2=T$.
The phases of $m$ and $c$ remain confounded; neither mean shift is separately
identified. General non-Gaussian increments identify $|m(q)|$, not universal
variance. Source-dependent changes, unstable within-window kernels,
code-dependent selection and correlated errors remain uncorrected.

Ratios of group means precede the product and modulus. Averaging groupwise
magnitudes targets something else and sampling bias is not zero. Four
Cartesian Fieller sets at marginal level $1-(1-L)/4$ give a Bonferroni
rectangle for joint level $L$, retaining the independent-group approximation.
Propagate all Cartesian interval components through the product/modulus;
retain weak-denominator, unbounded and unusable outcomes, and never clip
legitimate signed components. Conservatism can lose resolution. The real-valued
variant coverage study does not validate this complex extension, and real-only
retained moments cannot reconstruct discarded imaginary observations.
