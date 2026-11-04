# Four-panel kernel-robust coherence and its limits

## Population result

Retain independent, exogenous, permanent four-source codes and the existing
fourth-order source-diagonal extractor H4. Panels a,b belong to window A and
c,d to B. Conditional on counts, source nuisances and persistent positions,
disjoint-panel errors factor. Within each window, source s has kernel psi_sA or
psi_sB, arbitrary complex characteristic functions at q. They can change
arbitrarily between windows. Conditional on the complete kernel/count profile,
the displacement characteristic m=E exp[iq(D_sA-D_sB)] is common to all sources.
This last conditional requirement excludes motion/kernel correlation masquerading
as a common unconditional motion law. No Gaussianity is needed so far.

Write w_s=n_sa*n_sb*n_sc*n_sd (or its conditional/integrated version with the
same joint conditioning). Define

    RA = E Re H4(Za,conj(Zb),Nc,Nd)
    RB = E Re H4(Na,Nb,Zc,conj(Zd))
    C  = E [H4(Za,Nb,conj(Zc),Nd)+H4(Za,Nb,Nc,conj(Zd))
           +H4(Na,Zb,conj(Zc),Nd)+H4(Na,Zb,Nc,conj(Zd))]/4.

Then source-diagonal extraction and conditional factorization give

    RA=sum_s w_s |psi_sA|^2,
    RB=sum_s w_s |psi_sB|^2,
    C=m sum_s w_s psi_sA conj(psi_sB).

With positive RA,RB, Cauchy–Schwarz therefore proves

    G=|C|/sqrt(RA*RB)=|m| rho <= |m|,
    rho=|<sqrt(w)psi_A,sqrt(w)psi_B>|/(||sqrt(w)psi_A|| ||sqrt(w)psi_B||).

Equality holds exactly when the two weighted kernel vectors are complex
proportional on their contributing support. This includes arbitrary nonzero
source-uniform complex change. Source-specific changes can lower G; they need
not be small, centered, Gaussian, positive-real or calibrated individually.
The proof is an application of classical coherence/Gram-matrix geometry, not
a new inequality or general source-isolation principle.

If the actual displacement law additionally gives |m(q)|=exp(-q^2*T), then
T<=-log(G)/q^2. A lower confidence bound on G gives an upper bound on T.
It does not measure mean displacement: deterministic common inter-window shift
has |m|=1 even when images misregister. It also does not identify the variance of
an acquisition-wide random shift from many internal groups in one acquisition.

## What is not obtained

The bound is sharp as a statement about these moments. With two equal-weight
sources, psi_A=(1,1), psi_B=(exp(i phi),exp(-i phi)), rho=cos(phi) for
0<=phi<=pi/2. Any G<=|m|<=1 can be represented. Low G does not establish
nontrivial physical motion or an equivalence certificate. Population references
can vanish; finite estimates can be negative, exceed unit coherence or be weak.
None is clipped into the physical region to manufacture a result.

If physical increments are canceled by correlated imaging offsets, observed
coherence can equal1 despite nontrivial physical motion. Such correlation violates
the conditional factorization/common-law premise and is not repaired here.
With heterogeneous m_s, the same argument bounds G by max_s |m_s|; it does not
certify a source average or every source. Count-dependent source selection and
label-dependent contributions remain separate failures. The four-panel weights
and support differ from the old cubic weights under heterogeneous motion.
Sources missing a panel have zero four-panel weight, not recovered information.
All these limitations must accompany the output.

## Why this cannot remove arbitrary shared physical bias

For a recorded position Y_s(t)=X_s(t)+B_s(t)+epsilon_s(t), the transformation

    X'_s(t)=X_s(t)+U_s(t),  B'_s(t)=B_s(t)-U_s(t)

leaves every observed position, code, calibration read and count unchanged.
U can carry nontrivial physical motion. Thus the same complete observation law
can have different physical motion characteristics if arbitrary dependent
imaging bias is admitted. This is an observation-level ambiguity, stronger than
mere equality of the four moments. No function of those unchanged observations
can universally distinguish the models. An applicable independent constraint
on B, or an explicitly justified exclusion of this dependence, is needed.
The coherence theorem excludes such compensation through its conditional
model; the theorem is not empirical evidence that compensation is absent.

## A direct threshold witness

For g in[0,1], fixed phase theta and positive t, define group statistic

    W_g=Re(exp(-i theta)*C_g) - g*(t*RA_g+RB_g/t)/2.

If |m|<=g, then each independent group's expected W_g<=0: complex projection
is at most |C|, the coherence theorem bounds |C| by g sqrt(RA RB), and the
arithmetic–geometric mean bounds sqrt(RA RB) by (t RA+RB/t)/2. This avoids
discarding numerator/reference covariance by separately bounding four moments.
Theta and t may be fixed in advance or mixed over a prespecified finite set;
they cannot be selected on scored observations without valid adjustment.
Calibration on the same persistent source states is not automatically independent
training for a physical sampling claim. Calibration-only computations below are
therefore labeled design feasibility, never validation confidence statements.

For B_g=product_t total_count_gt, all four real component magnitudes are bounded
by B_g because the H4 code coefficients have absolute value<=1. Thus

    |W_g| <= B_g [1+g*(t+1/t)/2] = L_g.

Under independent groups with the per-group null mean condition, each
1+lambda*W_g/L_g is nonnegative and has expectation<=1 for0<=lambda<=1.
Products across groups, and a fixed mixture of such products over theta,t,lambda,
are e-values. Reject |m|<=g only at e>=1/alpha. This is a finite-sample bound
under the stated independent-unit model, not proof of that model. A t-test on
the same linear witness is a separate approximate diagnostic; successful normal
simulation coverage does not make it a distribution-free guarantee.

A simultaneous component Hoeffding rectangle is also implemented as a simple
outer check. Its common scaling preserves global count weights. Rectangles may
be much less informative than a covariance-aware witness. Neither finite method
can erase the physical interpretation gap by producing a smaller fitted error.

## Conditional code covariance

For code tuples of length r<=4, let h_r be1 for one unique code, -1/3 for two,
1/3 for three and -1 for four. The existing H4's code coefficients have exactly
these values. A quartet with j distinct physical source indices reduces to
h_j of those source codes, independent of the repeat positions. Each h_j is
canonical: averaging over any one of its independent source codes gives zero.
Different source subsets are therefore orthogonal. Exact enumeration gives

    E h2^2=1/3, E h3^2=1/6, E h4^2=5/24.

Let F(S) be the vector of collapsed quartic components if every source in S
shares one code, and let b_S=sum_(U subset S) (-1)^(|S|-|U|) F(U).
The diagonal mean is sum_i F({i}); for2<=|S|<=4 the b_S multiply h_|S|.
Consequently the exact conditional covariance is sum_S E[h_|S|^2] b_S b_S^T.
Real and imaginary components are stacked, retaining their dependence. This
extends the earlier source-subset calculation to this quartic statistic;
canonical decompositions and their polynomial coefficients are established.
Physical observation variation is an additional term, not this covariance.

## Strong information-rich alternatives

If a collision-free flag is known from the design for a group of m<=4 sources,
it has probability p=(4)_m/4^m. On that event every occupied code contains at
most one source, so sum_code F(code) is the source diagonal. Multiplying by
I(collision-free)/p gives the same expected moments. Colliding groups are
reported with zero contribution, not removed to select a favorable seed. This
Horvitz–Thompson construction can have very different variance from H4. It
requires independent design knowledge of the flag. Counting detected channels
is unsafe when sources may disappear; calibration-known occupancy alone does
not make that flag observable. No universal source-blind superiority follows.

An exact count-diagonal control variate D-W also has zero conditional code mean
when the true current four-count mass W is independently known. Replacing moments
by M-beta*(D-W) preserves their means for code-independent coefficients, but
W cannot be replaced by an unqualified calibration prediction. Calibration-only
variance reductions from such a control are not deployable inference results.
