# Conditional covariance of the count-matched Fourier contrast

## 1. Randomization and source-index decomposition

Condition on every finite source/panel phase sum and weighted count. Codes
c_i are iid uniform on {1,...,K}, K>=3, permanent across panels and independent
of these conditioned observations. For three code labels define h=1 if all
are equal, -1/(K-1) if exactly two are equal, and
2/((K-1)(K-2)) if all are distinct. Define
u(a,b)=(K*1[a=b]-1)/(K-1).

For any real or complex source coefficient tensor a_ijk, the aggregate statistic
is T_a = sum_ijk a_ijk h(c_i,c_j,c_k). Its exact decomposition is

    T_a = m_a + sum_{i<j} B_a,ij u(c_i,c_j)
                  + sum_{i<j<k} D_a,ijk h(c_i,c_j,c_k),
    m_a = sum_i a_iii,
    B_a,ij = a_iij+a_iji+a_jii+a_ijj+a_jij+a_jji,
    D_a,ijk = sum_{six permutations of i,j,k} a_permutation.

Proof: split the ordered source triples by support size. For support size two,
h(c_i,c_i,c_j)=u(c_i,c_j), regardless of which source is repeated. No identical
source contributions, Gaussian kernels, equal counts or panel independence are
needed for this *conditional randomization* statement.

E[u|one code]=0 and E[h|any two codes]=0. Thus distinct subsets, and the pair
and triple orders, are orthogonal, even if subsets overlap. Direct summation gives
E[u^2]=1/(K-1), E[h^2]=1/((K-1)(K-2)). Consequently

    Cov(T_a,T_b | observations)
      = sum_pairs B_a B_b / (K-1)
        + sum_triples D_a D_b / ((K-1)(K-2)).                 (1)

For real outputs this is ordinary covariance. For complex outputs, replace
B_b,D_b by conjugates for Hermitian covariance; also retain the unconjugated
pseudo-covariance, or apply (1) to the stacked real and imaginary outputs.
The latter avoids assuming circular complex noise. Equation (1) automatically
includes all cross-window, cross-angle and cross-frequency terms when the same
source codes are used. Averaging marginal variances would omit these terms.

## 2. Actual packaged statistic

Let z_it(q,angle)=sum_{events of source i,panel t} weight*exp(i*q*x_projected),
and n_it=sum weights. For each panel triple (p,r,s), define

    R: Re[z_ip * conjugate(z_jr) * n_ks],
    O1: Re[z_ip * n_jr * conjugate(z_ks)],
    O2: Re[n_ip * z_jr * conjugate(z_ks)].

Apply h to every ordered source triple. In each window the package averages
the two opposite-window count placements with weight 1/2. Its numerator is
(O1+O2)/2. Window A uses (0,1,2),(0,1,3); B uses (2,3,0),(2,3,1).
The arithmetic numerator/reference further average the two windows. Angle
averaging occurs *before* covariance, retaining code dependence across angles.
For K=4, the package's six-permutation Hadamard contraction has exactly the h
coefficients above. The general-K formula describes the corresponding symmetric
kernel; it does not alter the shipped four-code implementation.

Expand Re[x*conjugate(y)] as Re(x)Re(y)+Im(x)Im(y), or the analogous imaginary
expression. Each real component is a finite sum of separable real tensors
w_l x_li y_lj z_lk. Thus (1) applies directly to the implemented contributions,
including unequal counts and finite localization/event noise conditioned upon.

## 3. Physical variation and ratios

With physical observations Z independent of assigned codes,

    Cov(T) = Cov_Z(m(Z)) + E_Z[Cov_codes(T | Z)].              (2)

Additional alphabet size affects the second term, not the first under this
stipulated fixed-observation law. Shared acquisition motion is not independently
resampled by internal groups. Equation (2) needs the *actual* sampling unit.
Code-dependent offsets violate the conditioning and are outside (1)-(2).

For independent identically distributed groups and positive/nonzero population
reference, arithmetic a = mu_C/mu_R has first-order variance

    Var(a_hat) approximately Var(C-aR)/(G*mu_R^2).

This is the covariance *of the ratio influence*, not Var(R) alone. For the
geometric a=sqrt((mu_CA/mu_RA)*(mu_CB/mu_RB)), positive marginal ratios are
required. Its gradient is a/2*(1/mu_CA,-1/mu_RA,1/mu_CB,-1/mu_RB).
For log geometric attenuation omit the leading a; for log arithmetic use
(1/mu_C,-1/mu_R). For T=-log(a)/q^2 multiply the log gradient by -1/q^2.
These derivatives do not provide finite-sample ratio moments or guaranteed
confidence coverage. Keep historical Fieller sets and unresolved outcomes.
Different window weights can give arithmetic/geometric different targets.

## 4. Restricted uniqueness, not an information bound

Within fixed-coefficient trilinear aggregate estimators that are unbiased for
arbitrary fixed contributions, average over all code-name permutations. Code
assignment is invariant under that group; Jensen's inequality shows this cannot
increase scalar conditional squared error (and the matrix covariance difference
is PSD by applying the argument to all real linear combinations). Code-name invariance alone leaves FIVE equality-pattern coefficients:
alpha; beta_12,beta_13,beta_23 for the three possible repeated positions; and
gamma. It does not itself imply symmetry in the three panel arguments.
Arbitrary-contribution unbiasedness for one source gives alpha=1. For EACH
repeated-position pattern of two physical sources:
alpha/K+(K-1)*beta_position/K=0. Hence all three betas equal -1/(K-1).
For three independent source codes:
alpha/K^2+(K-1)*(beta_12+beta_13+beta_23)/K^2
 +(K-1)*(K-2)*gamma/K^2=0. These force the h coefficients. Thus h is the covariance-minimizing symmetric
member of this fixed-coefficient trilinear unbiased class, and symmetrization
maps every admissible member to h. It need not be uniquely minimizing before
symmetrization. This is a classical averaging argument, not universal optimality.
Nonlinear estimators, code-count dependent coefficients, biased estimators,
source identities and stronger physical assumptions are outside the class.

For count-only unit sources, B=6 and D=6. Therefore
Var(W)=18N(N-1)/(K-1)+6N(N-1)(N-2)/((K-1)(K-2)).
N=sum n_c is already observed in this toy, so W is not efficient for estimating
N without additional restrictions. For N=3,K=4, W is -1,3 or27, never zero.
If C=W*D, averaging C/W returns the ordinary mean of D, whereas the ratio of
sums can have variance inflation 17/3. This disproves an unrestricted coding
penalty, not the conditional covariance or restricted trilinear result.

## 5. Resource implications

For any fixed scalar influence v, write A=v'Cov_Z(m)v,
P=E sum_pairs (v'B)^2, Q=E sum_triples (v'D)^2. Then per-group
variance is A+P/(K-1)+Q/((K-1)(K-2)). At equal groups it decreases with K unless
P=Q=0. At B=K*G the relevant cost factor is

    f(K)=K*[A+P/(K-1)+Q/((K-1)(K-2))].                     (3)

Increasing K to L improves this stipulated budget exactly when f(L)<f(K).
The physical floor contributes K*A and can reverse a coding-only gain. For
large K and P,Q fixed, f(K) is approximately K*A+P+(P+Q)/K;
the toy sqrt(6N) rule is not a general optimum for actual Fourier observations.
If C=aR source-by-source and for every code assignment, the ratio influence
vanishes and there is no coding gain. Hypothetical alphabet changes must retain
counts, kernels, target and independent groups for this comparison; actual
hardware can violate each condition. Group splitting is not free.

## 6. Evaluation status

Derivation above precedes comparative outcomes. Numerical, retained-evidence
and literature results are recorded separately. This extension cannot validate
physical labeling exogeneity, stable kernels, independent source errors, the
unperformed randomized MERFISH intervention or independent molecular truth.

### Proof review correction (before manuscript promotion)

The first draft conflated code-name invariance with argument symmetry. The
five-orbit argument above repairs that omission: the three two-source
unbiasedness constraints, not invariance alone, force the repeated-position
coefficients to agree. The result and covariance implementation are unchanged.


## 7. Exact algebraic relation to third-order k-statistics

For code aggregates X_c,Y_c,Z_c and bar denoting the average over K bins,
our contraction is exactly K times

    k_111 = K^2/((K-1)(K-2)) *
       [bar(XYZ)-bar(XY)bar(Z)-bar(XZ)bar(Y)-bar(YZ)bar(X)
        +2bar(X)bar(Y)bar(Z)].

This is the classical unbiased multivariate third-cumulant polynomial (see
Schefczik and Haegele, Eq.29, linked in PRIOR_WORK.md). The identity is
algebraic; code-bin aggregates are dependent under fixed source assignment,
so the iid-bin cumulant sampling model is NOT our randomization law. The
source-index proof establishes the intended diagonal target under the actual
assignment law. The polynomial coefficients and index-partition technique
must not be claimed as new. The application-specific result here is the joint
source-conditioned covariance and its explicitly conditional design implications.

### Explicit four-to-eight-code decision boundary

For the SAME scalar target/influence and B=K*G, comparing (3) at K=4 and K=8
simplifies exactly to

    eight codes reduce first-order variance iff A < P/21 + 5Q/42.

Equality is a tie; the opposite inequality favors four codes. This boundary
includes the physical floor and numerator/reference cancellation. It is NOT
an observable source-blind recommendation, because P,Q generally require
unresolved source contributions or an independently justified source model.
The threshold's measurement meaning is application-specific; comparing two
known variance/cost formulae is elementary algebra, not a new design principle.


### Phase-preserving and orientation outputs

The same covariance applies to Re C_A, Im C_A, R_A, Re C_B, Im C_B, R_B;
references are real in the shipped phase-preserving method. For
log sqrt(|r_A*r_B|), the window-A gradient entries are
(1/2)*(Re C_A/|C_A|^2, Im C_A/|C_A|^2, -1/R_A), with the analogous B entries.
Zero complex numerators and nonpositive references remain unsupported.
`phase_log_gradient` is finite-difference tested and the underlying components
are checked against `phase_components` for every three-source four-code
assignment. `include_orientations=True` also exports both cross orientations;
their full covariance is verified against the historical `components` routine.
For log D, the gradient is (1/C_A,-1/R_A,-1/C_B,1/R_B). None of these derivatives
replaces the conservative shipped confidence-set construction or certifies
model assumptions from a quiet diagnostic.
