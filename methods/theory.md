# Supplementary Note 1. Calibration observables and their identifiability

This is the consolidated mathematical statement for the manuscript. It
supersedes conflicting statements in the chronological theorem/addendum
file, which remains an audit record. Algebraic identities, model-based
identification and numerical classifications are distinguished below.
The Hadamard identity, nuisance-space projection and affine registration
are established mathematical tools. Fourier dynamical normalization and
temporal source-diagonal cumulants also have established microscopy precedents
(main-text references [17–19]). The construction here is the conditional
matched-count observable under permanent source-code randomization, not a new
general cancellation principle or a demonstrated transport/reconstruction advance.

## 1. Random-code source-diagonal extraction

Let sources s have fixed positions θ_s and independent uniform codes
c_s in {1,2,3,4}. Let h_1,h_2,h_3 be the nonconstant rows of the four-point
Hadamard matrix, with signs chosen so that h_1 h_2 h_3 = 1 pointwise.
For source contributions a_s,b_s,d_s, define A_j = Σ_s h_j(c_s)a_s and
analogously B_j,D_j, and

    H₃(a,b,d) = (1/6) Σ_{permutations (j,k,l) of (1,2,3)} A_j B_k D_l.

**Lemma 1.** Conditional on all contributions, if the codes are independent
of them, E_c H₃ = Σ_s a_s b_s d_s.

**Proof.** Expand each product over source indices (s,t,u). When all three
indices differ, independence gives a product of zero-mean signs. When two
coincide, the remaining source contributes a zero-mean sign. When all
coincide, the product of the three distinct Hadamard rows is one. Each
of the six permutations therefore has expectation Σ_s a_s b_s d_s. □

This is an expectation over physical source-code randomization. Randomly
splitting events after acquisition is not an implementation: events from
one persistent source must carry the same code. A realized finite field
has randomization variance; this lemma does not say its estimate is exact.

## 2. Matched-count calibration identity

Project positions onto a unit direction and use spatial frequency q.
Within window w, source s has persistent displacement D_sw and event
error ε_stℓ. An event in panel t of window w(t) is

    X_stℓ = θ_s + D_s,w(t) + ε_stℓ.

Write Z_st(q)=Σ_ℓ exp(iqX_stℓ) and n_st for its count. The aggregate
code moments are formed from these source contributions without observing
source identity. Use disjoint event panels a,b in window A and c in B.

The following assumptions are sufficient:

1. Source codes obey Lemma 1 and are independent of positions, counts,
   persistent displacement and event errors.
2. Conditional on counts, source nuisance and persistent displacements,
   errors in disjoint panels factor, with the same source-specific
   characteristic function ψ_s(q), independent of the displacement
   realizations. This allows heterogeneity between sources but excludes
   panel-dependent kernels, shared events and displacement-dependent
   event-error kernels.
3. Conditional on the quantities determining the count-product weights
   and ψ_s, E exp[iq(D_sA−D_sB)] = a(q), the same real positive attenuation
   for every source. Count-dependent displacement laws require a different
   estimand and are not covered merely because the counts are balanced.
4. The common denominator below is nonzero. Finite-sample ratios must also
   satisfy the estimator's sign and conditioning checks.

Define population moments

    R_A = E Re H₃(Z_a, conjugate(Z_b), n_c),
    C_A = E Re {H₃(Z_a,n_b,conjugate(Z_c))
                   + H₃(n_a,Z_b,conjugate(Z_c))}/2.

**Theorem 2.** R_A = Σ_s E[n_sa n_sb n_sc |ψ_s(q)|²] and
C_A = a(q) R_A. For fixed source kernels, write w_s=E[n_sa n_sb n_sc]
to obtain R_A=Σ_s w_s |ψ_s(q)|². Thus C_A/R_A=a(q).

**Proof.** Lemma 1 removes every cross-source term. In the reference,
θ_s and D_sA cancel. Conditional panel independence gives |ψ_s|².
In either cross-window term θ_s cancels, the error factor is again
|ψ_s|², and Assumption 3 contributes a(q). All three terms use exactly
the same product of three counts, including the count in the unphased
panel. Summing over sources gives the claim. □

**Why matching the count product matters.** The issue is not removed merely
by dividing each moment by its total number of source pairs. To see the
distinction, condition on fixed source counts and stable kernels at a nonzero
frequency, and retain the common attenuation assumption above. Write
K_s=|ψ_s(q)|². A source-diagonal within-window pair moment has weights
n_sa n_sb; a cross-window pair moment has weights n_sa n_sc. Suppose both
total pair counts are positive, and normalize each set of weights to sum to one:

    u_s = n_sa n_sb / Σ_j n_ja n_jb,
    v_s = n_sa n_sc / Σ_j n_ja n_jc.

The corresponding normalized population moments obey

    R₂ = Σ_s u_s K_s,     C₂ = a(q) Σ_s v_s K_s,
    C₂/R₂ = a(q) [Σ_s v_s K_s / Σ_s u_s K_s].

The bracket need not equal one, even when every source kernel is stable
between windows and there is no motion. For this particular ratio to equal
a(q) for every positive source-kernel vector, u=v is necessary and sufficient.
Sufficiency follows by cancellation. For necessity, equality for every K
implies Σ_s(v_s−u_s)K_s=0; varying one positive coordinate at a time forces
v_s=u_s for each source. Admissible independent Gaussian event errors already
allow such variations within 0<K_s<1 at a fixed nonzero q.

Theorem 2 puts the same triple count product in every term, including the
otherwise unphased panel. This removes the weighting factor without learning
each K_s, equalizing the counts or assigning sources. It identifies the common
attenuation on the supported triple-product population; it does not measure
sources absent from that support. This algebra explains the construction,
not a universal requirement for third-order measurements or random coding.
Homogeneous or known kernels, matched pair weights, or another observation
model can permit other valid normalizations. In particular, the calculation
is not an empirical failure test of kICS, DDM or SOFI. Panel-dependent kernels,
code-dependent selection and heterogeneous motion remain outside Theorem 2.
Nor does this factor establish the cause of all empirical reference gaps:
those comparisons also change source selection, weights and temporal support.

The heterogeneous-sample kICS extension [34] already handles an immobile blinking population and removes shared image-PSF dependence through a normalized spectral autocorrelation (its Equation 14), with a stated transport/photophysical model (Equation 15). This precedes any claim here about accommodating heterogeneity. The count identity instead matches panel-dependent count products against source-specific localization-error kernels under permanent source codes. It does not establish that kICS fails on its own target, or that random codes are necessary for every dynamical measurement. No head-to-head performance comparison was run.

Swapping A and B gives the second matched ratio, with potentially
different source weights but the same a(q). Arbitrary count heterogeneity
is allowed subject to the conditional factorization above; unbalanced
second-order ratios generally use different weights and do not enjoy this
identity. A ratio of sample moments need not be unbiased at finite sample
size. Group/source resampling is needed for uncertainty, with the unit
matching the dependence structure. These are population identities. If a
displacement component is shared by all groups in one acquisition, internal
group resampling conditions on that realization and cannot quantify its
across-acquisition variation. More groups do not create independent
realizations of a shared shift or identify its ensemble variance from one
window pair. Independent realizations or a specified shared-motion model are
needed for that additional interpretation.

If D_sw=G_w+U_sw, the two windows are independent, and the projected
components G_w and U_sw are independent zero-mean Gaussians with per-window
variances γ² and τ²,
then Var(D_sA−D_sB)=2(γ²+τ²), so

    a(q)=exp[−q²T],   T=γ²+τ²,   T=−log(a(q))/q².

Independence and a finite variance alone do not justify the exponential.
Without Gaussianity the contrast identifies the displacement attenuation
at the measured q values. Recovering variance from curvature at q=0
requires the appropriate second moment and a justified limiting estimator.

## 3. What the split diagnostic does and does not prove

Let r_A=C_A/R_A and r_B=C_B/R_B. Under Theorem 2, D=r_A/r_B=1.
Rejection is evidence against the model; non-rejection is not a certificate
of its validity. Different violations can cancel in D. The repository
contains a population witness with D=1 and frequency-dependent bias
(`outputs/d-limit-sep18`), and a >40 nm² stability-violation witness.
Frequency agreement and orientation checks are additional necessary
checks, not a complete family of sufficient tests.

For the special between-window change ψ_sB(q)=c(q)ψ_sA(q), identical
across sources, the real-valued contrasts give

    r_A = Re(c) a,     r_B = Re(c) a / |c|².

When c is positive real, sqrt(r_A r_B)=a, although D can reject.
For arbitrary complex c, the geometric mean is a|Re(c)|/|c| and generally
fails to recover a; a mean shift is a counterexample. Positivity of the
ratios is required for the logarithmic estimator. The positive-real
qualification is essential, not an implementation detail.

### 3.1 Explicit temporal support and a qualified fourth-order alternative

The four panels do not measure one instantaneous lag. Outside the constant-window
model of Section 2, define the temporal target before interpreting a fitted
Gaussian parameter. For a source with panel counts n_it, let W_i=prod_t n_it
and mu_it be its mean retained event coordinate in panel t. A directly defined
signed target is

    B_obs = sum_i W_i (mu_i1-mu_i4)(mu_i2-mu_i3) / (2 sum_i W_i).

This weights only sources observed in all four panels. In the observation model
X_it=theta_i+D_it+epsilon_it, conditionally zero-mean errors with zero cross-panel
covariance cancel from its expectation. Their source/panel variances may differ.
The latent target is then half the weighted expectation of
(D_i1-D_i4)(D_i2-D_i3). It is half a covariance only when the increments are
centered; otherwise it also contains their mean product. Event timing, selection
and count conditioning are part of the assumptions. Finite panels average actual
event-time pairs; replacing them by midpoints changes the target in general.

The already-developed fourth-order code contrast supplies this target in
expectation. One explicit H4 kernel sums products of four per-code contributions,
with coefficients 1, -1/3, 1/3 and -1 when their code tuple has respectively
one, two, three or four distinct labels. Under independent uniform four-code
labeling, averaging this kernel gives one for a single source index and zero
for each other source partition. This follows by enumerating the possible
coincidences of the independently drawn labels. Thus H4 extracts the fourth
source diagonal without assignments. For code counts N_t and coordinate sums S_t,

    D0 = H4(N1,N2,N3,N4),
    N0 = [H4(S1,S2,N3,N4)+H4(N1,N2,S3,S4)
          -H4(S1,N2,S3,N4)-H4(N1,S2,N3,S4)]/2.

The ratio of their code expectations is B_obs; the realized ratio is not
unbiased. The bounded alternative uses phase sums Z_t(q) and

    C4(q)=E Re H4(Z1,Z2,conjugate(Z3),conjugate(Z4)),
    R4(q)=E Re H4(Z1,conjugate(Z2),Z3,conjugate(Z4)).

For u=X1-X4 and v=X2-X3, the phase patterns are u+v and u-v. Their squared
difference is 4uv, proving, with finite second moments and nonzero count mass,

    lim_(q->0) -log[C4(q)/R4(q)]/(4q^2) = B_obs

for population moments of the retained observations. Finite-frequency estimates
need additional assumptions. Independent symmetric event errors contribute the
same real product of four characteristic factors to both patterns, even if each
source/panel has a different kernel. With a centered jointly Gaussian latent
four-coordinate law and a common conditional B across sources, this yields
C4/R4=exp(-4q^2 B). This qualified robustness does not follow for arbitrary
asymmetric kernels, changing means or cross-panel dependence. It is not a
general repair of nonstationarity or an efficiency guarantee.

For a common centered process with covariance K and panel-time sampling
densities p_t independent of source, B=(Kbar_12+Kbar_34-Kbar_13-Kbar_24)/2,
where Kbar_ab integrates K(s,t) against p_a(s)p_b(t). Independent constant
windows recover T in Section 2. Brownian motion with increment variance
2 kappa h gives B=kappa(E t3-E t2) for ordered nonoverlapping panels. Linear
deterministic drift can give positive B despite zero stochastic variance;
oscillatory covariance can give negative B. Neither interpretation is a universal
resolution or molecular-distance error bar.

An unrestricted panel-dependent error mean can absorb the displacement without
changing any observation. Correlated errors can contribute to the same temporal
product. These ambiguities require additional assumptions or independent
information, not more precise estimation of the same moments. Denominator
uncertainty also matters: unbounded or disconnected Fieller ratio sets must not
be replaced with finite intervals by clipping estimates or omitting invalid
resamples. The development derivation, exact checks and retained negative results
are in `outputs/temporal-definition-sep28`; this alternative is not promoted
to the production estimator or to independently validated molecular calibration.

### 3.2 Phase-preserving extension for a source-uniform complex kernel change

The positive-real restriction of the historical geometric construction can
be relaxed by retaining the complex third-order cross moments. This extends
the specified count-matched observable; complex Fourier phase and magnitude
are established tools (main-text reference [17]). It is not a new general
coherence principle or an empirical superiority claim.

Keep Theorem 2's code exogeneity and conditional panel factorization. Kernels
must be stable within each window, but let ψ_sB(q)=c(q)ψ_sA(q), with the same
nonzero complex c for every included source and group. Let
m(q)=E exp[iq(D_sA−D_sB)] be common conditional on the count/kernel weights;
it need not be real. Require both window references to be nonzero. Define
complex cross moments C̃_A and C̃_B using the same H3 phase/count placements
as before, without taking their real parts. The within-window references
remain real. Their source weights can differ between windows.

**Corollary 2.1.** Under these assumptions,

    C̃_A/R_A = m conjugate(c),
    C̃_B/R_B = conjugate(m)/conjugate(c),
    (C̃_A/R_A)(C̃_B/R_B) = |m|².

Consequently sqrt(|(C̃_A/R_A)(C̃_B/R_B)|)=|m| independently of c.

**Proof.** Source-diagonal extraction is complex multilinear. Within window A,
the error factor is |ψ_sA|². Its cross-window counterpart is
ψ_sA conjugate(ψ_sB)=conjugate(c)|ψ_sA|² and its displacement factor is m.
Both terms have the same triple count product, so the first ratio follows.
In window B the cross error factor is c|ψ_sA|², the displacement factor is
conjugate(m), and the reference factor is |c|²|ψ_sA|². Dividing gives the
second ratio regardless of its different source weights. Their product
cancels c. □

This is an information-preserving extension. At one frequency, the retained
5-nm mean-shift witness with T=9 has the same real reference/cross moments as
a stationary-kernel Gaussian model with T=22.058... at q=.10. Thus those
real moments alone do not distinguish the two models. Their complex moments
do. Keeping phase recovers 9 in the original population witness at both
retained frequencies, rather than the historical 22.058/22.884. The original
failure remains unchanged. This is a closed-form and exact-enumeration
control, not another empirical acquisition or a coverage study (Source Data 43).

For a Gaussian displacement increment with mean µ and variance 2T,
m=exp(iqµ−q²T), so −log|m|/q²=T without requiring µ=0. The phase of m is
confounded with the phase of c: neither mean displacement nor kernel shift
is identified separately. Internal groups condition on a shared acquisition
realization and cannot identify its across-acquisition variance. General
non-Gaussian increments identify |m(q)|, not a universal variance or molecular
distance error. Source-dependent changes, within-window instability,
code-dependent selection and correlated errors are not repaired.

The sample estimator takes ratios of group means and then the product and
modulus; averaging groupwise magnitudes would target a different quantity.
Its sampling bias is not zero. Four component Fieller sets, each at
1−(1−L)/4 for joint confidence level L, give a Bonferroni rectangle for the two
complex ratios. No independence between components/windows is required;
marginal coverage still relies on the independent-group approximation.
Propagate every Cartesian interval component through the modulus and product,
retaining weak-denominator, unbounded and unusable outcomes. Signed Cartesian
components are legitimate and must not be clipped. Extra conservatism can
lose useful resolution; the Source Data 71 coverage results concern the
historical real-valued variants and do not validate this extension's sampling
performance. Real-only retained moment tables cannot reconstruct discarded
imaginary observations, and no new cohort score or automatic selector is supplied.

## 4. Fixed-code mixture nuisance and the correct weighted null space

Let B be a fixed binary codebook with G codewords. On its E covered
bit-pair edges, N_eg is one when codeword g covers edge e. Let κ_g>0
be abundance/detection weights and D=diag(Nκ). In the local, normalized
mixture model the nuisance Jacobian is

    L = D⁻¹ N diag(κ).

For example, linearizing −log[Σ_g N_eg κ_g exp(−q²σ_g²)]/q² at q=0
produces (Lσ²)_e, after removing the intercept. At nonzero q or a
different expansion point, the Jacobian and weighting must be specified
again. The statement below is local identifiability of this stated model,
not a global theorem about arbitrary images.

**Proposition 3.** rank(L)=rank(N) and

    ker(Lᵀ)=D ker(Nᵀ),    dim ker(Lᵀ)=E−rank(N).

**Proof.** Both diagonal factors are invertible on covered edges, proving
rank equality. Lᵀu=0 iff NᵀD⁻¹u=0, proving the null-space formula. □

For y=Hβ+Lη, a linear edge contrast uᵀy is free of η iff uᵀL=0.
A target component of β is identifiable only if it is recoverable from
the projected design P H, where P projects onto ker(Lᵀ); a nonzero
nuisance-free edge contrast alone does not prove all target parameters
identifiable. With additional nuisance Zα, project against [L Z].
Full row rank of N gives no nonzero mixture-free edge contrast. Rank
deficiency gives possible contrasts, whose directions depend on κ under
the normalized convention. Unknown κ therefore matters even when the
dimension is known. For unnormalized sums L=N diag(κ), the null space is
instead ker(Nᵀ).

Zero-abundance codewords violate the invertibility assumption. For a specified
active population, remove their columns and any uncovered edges before
applying the rank statement. A full-book classification need not describe
the identifiable contrasts of every restricted expressed-gene population.

Since rank(N)≤G, the null-space dimension is at least max(0,E−G).
Pair balance alone is neither a sufficient nor necessary rank condition.
For an integer incidence matrix, rank modulo a prime is a lower bound
on real rank: a nonzero modular minor has a nonzero integer determinant.
Modular full row rank therefore certifies real full row rank. A deficient
modular rank alone does not certify real rank deficiency. Floating-point
ranks and the historical 0.05 projection cutoff in the 33-codebook survey
are numerical diagnostics. Its 20/22 MERFISH-family
Class C count combines exact-rank and low-sensitivity cases and must not
be called 20 proofs of exact blindness. The 0.05 threshold is not a
mathematical boundary. Assignment-ensemble variance floors are empirical
simulation results for stated gene pools, not exact universal constants.

## 5. Plane-factorizable blur cannot be separated without more information

Suppose the observable log amplitude contains

    α_i(q)+α_j(q)−q²(T_i+T_j)/2 + mixture terms.

**Proposition 4.** With unrestricted plane functions α_i, node-additive
T_i is not identifiable: T_i→T_i+δ_i and
α_i(q)→α_i(q)+q²δ_i/2 leave all amplitudes unchanged.

This follows by substitution. An all-order version holds whenever the
admissible plane kernel can absorb the displacement characteristic
factor; positivity/physical admissibility must hold. The local argument
does not require asserting that every arbitrary transformation remains
a valid characteristic function. The fixed-code projection certificate
does not remove this additional confound by itself.

## 6. Randomized-round retrofit: conditional contrasts, not empirical sufficiency

For a constant-weight h>1 codebook, define K_i=Σ_g B_gi κ_g and
K_ij=Σ_g B_gi B_gj κ_g. An independent per-molecule randomized readout R
samples each bit's population with the same probability. Under stable
detection/mixture weights and the additive second-moment model, let μ_ij
be the pair's abundance-weighted nuisance mean and μ_iR the bit-i mean.
Then w_ij=K_ij/[(h−1)K_i] obeys Σ_j w_ij=1 and

    Σ_j w_ij μ_ij = μ_iR.

**Proof.** Each codeword containing i contributes its κ_g nuisance to
exactly h−1 pairs (i,j), cancelling the denominator. □

Consequently m_iR−Σ_j w_ij m_ij = T_R−Σ_j w_ij T_j, with T here
including any inseparable plane blur. Let R be the resulting linear map
from all bit budgets and T_R to these contrasts. R1=0; recovering every
relative budget additionally requires rank(R)=number of original bits.
That rank must be checked for the specified book/weights; it is not
guaranteed merely by adding a randomized channel. Estimated abundance
weights and finite randomization introduce sampling uncertainty.

A replicated plane contributes information only through a specified
replicate noise and motion model. One replicate does not automatically
identify every free plane blur parameter. For an augmented model y=Xβ+Zη,
an estimand Cβ is identified iff every (δβ,δη) satisfying
Xδβ+Zδη=0 also satisfies Cδβ=0. This elementary kernel criterion is the
appropriate sufficiency check. The numerical implementation `codedcal.identifiability.assess` checks this
criterion with an explicit SVD tolerance and returns either observation
weights or an invisible parameter perturbation. It checks a supplied model,
not the physical validity of that model. No empirically randomized MERFISH
acquisition exists here.
The shipped retrofit estimator and its noise/rank demonstrations are
simulation only; they do not establish the proposed acquisition's success.

### 6.1 Sufficient alternatives and an explicit failed design

Let b denote effective plane blur in budget units and suppose R has rank n−1
on n budgets (including R), with kernel span{1}. The contrast observation is
c=R(T+b). Suppose additional observations are v=H b, where the physical
meaning and independence of these observations must be justified. Let C be
the n−1 contrasts to the last plane.

**Proposition 5.** All contrasts CT are identifiable from (c,v) iff
ker(H) is contained in span{1}.

**Proof.** An invisible perturbation satisfies δT+δb=a1 and Hδb=0.
Consequently CδT=−Cδb. This vanishes for every invisible perturbation iff
every vector in ker(H) is constant. Conversely any nonconstant δb in
ker(H), paired with δT=−δb, is a witness to non-identification. □

Thus H=C suffices: relative blur information is enough; the absolute blur
level need not be known. One measurement H=e_p^T does not suffice for n>2.
For any unmeasured j, δb=e_j and δT=−e_j leave both observables unchanged
and alter a budget contrast. These perturbations can be arbitrarily small
around strictly positive budgets and blur, so non-negativity does not remove
the ambiguity. A repeated image is not automatically an observation of b;
correlated localization errors and within-repeat motion require a model.

Alternatively, if b=a1 is known a priori to be common, Rb=0 and budget
contrasts are identifiable without blur observations. This is a substantive
physical assumption. For a declared lower-dimensional blur basis B, replace
b by Bα and test the augmented matrix [R,RB;0,HB] against [C,0].

Constant codeword weight alone also does not guarantee rank(R). The
four weight-two words 1100, 0110, 0011, 1001 with equal abundance give a
rank-two map on five budgets. The software refuses a full-budget inverse in
this case. Computational recovery in a full-rank book and the sufficient
models above remain design-time results, not instrument validation.

### 6.2 Abundance-path rank loss and conditional precision

The rank certificate is pointwise in abundance, not uniform over populations.
Write S=B^T diag(κ) B, D=diag(B^T κ) and H=S−D. With positive bit
abundances, the weighted partner matrix is P=D^−1 H/(h−1), and P1=1.
Therefore 1 belongs to col(P), so rank([-P,1])=rank(P)=rank(H).
Along a convex path of nonnegative abundance vectors with positive bit
coverage, H varies continuously and remains symmetric. A change in its
number of negative eigenvalues between endpoints entails a zero eigenvalue
on the path, and hence loss of relative-budget information. Full-rank
endpoints do not rule this out. This is an elementary spectral consequence,
not a new rank theorem.

The exposed-data design audit locates such paths for both the original and
all twenty donor-optimized assignments (`outputs/root-causes-oct01`). It
explains why optimizing inverse noise at finitely many donor populations
cannot certify robust precision between them. `codedcal.precision.abundance_path`
locates sign-change failures; a negative screen does not certify a path or
convex hull. Numerical endpoint checks use the estimator's rank tolerance.
For a justified contrast covariance V, a specified relative estimand C has
conditional covariance C R^+ V (R^+)^T C^T. The opt-in precision check requires
an explicit SD tolerance; it does not infer omitted noise, bias or physical
randomization. All these retrofit calculations remain simulation-only.

## 7. Spatial displacement is a separate observable

Assume informative plane pairs provide unbiased local displacement
observations with a common, correctly specified spatial support model:

    d_ij(x̄)=τ_i−τ_j+(A_i−A_j)(x̄−c).

**Proposition 6.** If each required edge has a full-column-rank spatial
design [1,x̄−c] and the informative comparison graph is connected, the
relative affine coefficients are identifiable up to a common affine gauge.

**Proof.** Full spatial rank identifies the edge's intercept and gradient
difference by linear regression. A connected graph's incidence operator
has only constant node vectors in its kernel. Apply this separately to
each affine coefficient. □

In two spatial dimensions a general affine fit needs at least three
non-collinear informative locations. An isotropic scale fit has fewer
parameters but still requires a nonzero spatial lever arm. Disconnected
graphs have separate gauges. This proposition is conditional on valid
displacement observations: repeated textures, tissue gradients, weak shared
content or density weighting can invalidate an image correlation peak as
an unbiased displacement. The correlation estimator is not guaranteed by
the codebook alone. Code randomness is unnecessary for this conditional
registration result; standard affine-registration mathematics supplies
the proof. The demonstrated capability is estimating useful chromatic
scale before barcode decoding, with measured coverage and failure modes.

## 8. Joint interpretation through a specified measurement

The spatial and temporal results constrain different terms of measurement
error. Let e be the vector of coordinate errors after an applied spatial
correction, with mean m and covariance K under a declared source population,
panel schedule and sampling scheme. For a fixed linear measurement L,

    E ||L e||² = ||L m||² + tr(L K Lᵀ).

This standard bias–covariance identity provides the connection, not a new
estimator. Spatial displacement observations constrain parts of m. Temporal
contrasts constrain attenuation functions or, under additional assumptions,
second-moment functionals; they do not identify all of K. K includes sampling noise, cross-source dependence and uncertainty in
estimated correction maps, with cross terms unless independence is justified.
A nonlinear molecular distance or decoding probability requires an additional
measurement model. Fixed label-to-anchor offsets are not removed by the identity.

For two sources with four panels, put u_s=e_s1−e_s4 and v_s=e_s2−e_s3.
For centered errors let B_s=E[u_s v_s]/2 and
B_12=E[u_1 v_2+u_2 v_1]/4. Then

    B_1+B_2−2B_12 = E[(u_1−u_2)(v_1−v_2)]/2.

This is a signed relative-change band functional, not the variance of a final
separation. Nonzero temporal means contribute their products. The existing
fourth-order self/distinct code construction supplies related aggregate
observables under its randomization and acquisition assumptions. Its self and
pair weights generally differ; homogeneous conditional models or explicitly
matched weights are needed before combining them. Aggregate pair covariance
cannot be substituted for a selected neighboring pair's covariance. No new
finite-sample or empirical accuracy guarantee follows from this identity.

An exact source-level information audit uses eight panel coordinates and 36
free symmetric covariance entries. Granting two individual-source bands plus
the distinct-source band identifies the relative band above, but not the
variance of the difference of four-panel source averages. Even granting the
full covariance of six within-source increments leaves that average-separation
variance unidentified. Add independent constant source offsets of variance
kappa: every temporal increment is unchanged, while the average-separation
variance increases by 2 kappa. This is a positive-semidefinite ambiguity.
Independent stationary panel noise can also change the average's variance
without changing these particular temporal band functionals.

Thus successful spatial correction need not reduce a temporal contrast, and
a temporal contrast need not predict residual molecular accuracy. The combined
framework reports identifiable task components and unresolved ones; it does
not turn MERFISH spatial correction into empirical validation of randomized
temporal calibration. Exact rational rank checks, six covariance witnesses
and a retained-data feasibility audit are in
`outputs/joint-error-model-sep29` (Source Data 61). No complete independent
joint molecular-risk comparison is established by the retained artifacts.
