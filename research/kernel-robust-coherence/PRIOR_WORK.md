# Closest-work boundary

Targeted primary-source check, 29 October 2026. This is not exhaustive priority
clearance. No cited paper evaluated this implementation.

* Nieuwenhuizen et al., [Image Resolution in Optical Nanoscopy, 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4149789/),
  Eq.1 and the repeated-localization correction: normalized spectral correlation
  and the need for independent subdivisions precede this work. Our numerator
  and within-window references first use permanent codes to isolate physical
  sources with matched four-count weights. That is a different construction;
  neither normalized coherence nor Cauchy–Schwarz is a new contribution.
* Kolin, Ronis and Wiseman, [k-space Image Correlation Spectroscopy, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1578478/),
  Theory and Results: separating transport from fluorophore photophysics and
  optical factors is established, with experimental microsphere/cellular
  demonstrations. We do not have equivalent independent physical validation.
  The narrower extension here is a one-sided, count-matched coherence under
  arbitrary source-specific between-window kernel changes, conditional on
  within-window factorization and a common displacement characteristic.
* Sehayek et al., [Rapid ensemble measurement of protein diffusion and probe blinking dynamics in cells, 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC9680803/),
  especially Eqs.14–15: photophysics-aware normalized reciprocal-space inference
  handles heterogeneous immobile populations with its own model and data.
  We have not demonstrated performance superiority to that workflow.
* Culley et al., [Sub-diffraction error mapping for localisation microscopy images, 2021](https://www.nature.com/articles/s41467-021-25812-z):
  the publisher's indexed text explicitly describes shared reconstruction bias
  in split-data FRC. Full publisher text was inaccessible in this audit. This
  reinforces the existing shared-bias limitation; internal agreement is not an
  independent physical accuracy reference.
* Canonical polynomial covariance, unbiased cumulants, randomized moment
  sketches and uncertainty-guided design already appear in the preceding
  [finite-code closest-work audit](../finite-code-covariance/PRIOR_WORK.md).
  H4 is the existing extractor. The exact source-subset covariance here is an
  application-specific calculation, not a new general decomposition.

The collision-aware comparator uses ordinary inverse-probability weighting.
The direct threshold witness uses classical complex projection and AM–GM;
the finite test follows by multiplication of independent nonnegative factors
with null expectation at most one. Those elementary arguments are proved in
THEORY.md and tested directly rather than claimed as new statistical theory.

The substantive deliverable is therefore a narrowly qualified robustness target,
its exact coded-moment covariance, and a calibration-only design improvement.
None of these establishes the original physical equivalence claim, a novel
universal resolution metric, chemical intervention validity or a universal win
against existing correlation microscopy.
