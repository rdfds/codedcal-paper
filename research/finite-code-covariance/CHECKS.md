# Verification record and repairs

- Core numerical/analytic controls: 34 passed with Python 3.13 in the existing
  validation environment. Includes exact rational count moments, independent
  source-index reference, exhaustive general-K moments, all 64 Hadamard kernel
  entries, actual orientation and phase outputs, Re/Im covariance, permutation,
  scaling, gradients, total covariance and explicit model counterexamples.
- Shipped estimator: all 232 packaged tests passed (82.73 seconds). No shipped
  numerical or user-facing module was edited by this extension.
- Historical snapshot: reproduce/verify.py passed; all 19,089 imported files /
  5,254,683,652 bytes and original scientific evidence remain unchanged. See
  historical-integrity.log. This is integrity, not physical validation.
- Numerical evidence: verification.json independently recounts all planned
  cases. verification-all-intervals.json additionally recalculates every one of
  the 259,200 confidence sets independently from retained sufficient statistics.
  There are also 1,296 independent full-array API checks.
- Final deterministic replay: replay-equality.json compares every world file
  against the first successful covariance run, without selecting outcomes.

## Corrections retained explicitly

1. The initial optimality proof treated code-label invariance as argument
   symmetry. It actually leaves five equality orbits; separate two-source
   unbiasedness equations force the three repeated-position coefficients to
   agree. Fixed in the proof and tested before promotion. Covariance unchanged.
2. A new shared-motion control used an arbitrary >100-fold assertion, which
   was false for the chosen toy's noise level. Its exact total-covariance check
   passed. Replaced that unsupported threshold with the exact shared-versus-
   independent sampling identity and non-vanishing physical-floor checks.
   No estimator, scientific outcome or qualification threshold was changed.
3. The first report assembly was invoked before the finite worker completed
   and correctly refused to summarize. The empty output directory was removed;
   completion checks now precede directory creation. No partial report was
   presented as a completed result.
4. The first-pass resource scores mix a prefix coding estimate with a full-
   population floor. They are retained as preliminary descriptive estimates.
   The separately declared finite-budget check uses one consistent sixteen-
   state observation law for mean, both variance terms and each variant target.
5. Monte Carlo variance discrepancies remain explicit (23/864 exceed three
   estimated MC standard errors); extra draws were not added to improve them.
   Small-case exact identities and the general proof carry different evidence
   from a finite numerical covariance approximation.

No claim of external registration, independent physical confirmation, universal
optimality, exact Fieller coverage or general novelty follows from these checks.
