# Finite-code covariance and design extension — results

27 October 2026. Bounded conditional theory/computation; original manuscript,
Source Data 69–72, failed results and shipped estimator remain unchanged.

## Main finding

The exact covariance of the implemented source-coded Fourier moments separates
into pair and triple Gram matrices with alphabet factors 1/(K−1) and
1/((K−1)(K−2)). Physical-field variation is a separate floor. For a fixed
scalar ratio influence, eight codes beat four at B=K×G exactly when
**A < P/21 + 5Q/42**. These quantities require source contributions or a
justified source model; this is not an observable source-blind recommendation.

There is a real **model-conditional** variance/design improvement, but no
universal improvement in usable measurement. Fewer independent groups can
make Fieller sets unusable despite lower first-order variance. No physical
hardware cost, new random-readout experiment or independent molecular truth
has been established.

## Verification and complete scope

- All 126 historical worlds: 2,016 known-source and 4,032 coded prefix
  contribution vectors reproduced bit-for-bit.
- All nine views/eight acquisitions and both noise laws included.
- 108 exogenous worlds receive conditional covariance/design evaluation;
  all 18 code-dependent-offset worlds remain explicit theorem exclusions.
- 512 independent code assignments per fixed sixteen-group observation
  prefix check covariance; these are numerical replicates, not acquisitions.
- Median relative matrix discrepancy is 7.9%; maximum 28.6%.
- 23/864 component variance differences exceed three *estimated* Monte Carlo standard errors; maximum 5.38.
  These discrepancies are retained; estimated Monte Carlo errors are not
  treated as calibrated Gaussian agreement tests. No extra draws were added
  after disagreement. Exact proof/enumeration is a separate check.
- 64,800 finite-budget simulations (200 per K/world), K=4,8,16 and
  B=256, give 259,200 planned frequency/variant cases. No K is selected
  from the outcomes. Every sign/support failure and confidence-set type
  is retained. Original failures are not replaced by these simulations.

## Equal-budget results

The observation law is uniform on the SAME fixed sixteen reconstructed
groups, with fresh independent permanent codes. Both terms of total
covariance and each variant’s target use that same finite-state population.
Counts, photons and kernels are held unchanged by stipulation as K changes.
This establishes a conditional computation, not a measured acquisition design.

| Variant | K / groups | Median predicted variance / K4 | Conditions with lower variance | Median available-ratio RMSE / K4 | Usable + covering / all planned calls | Conditions better / worse / tied usable coverage |
|---|---:|---:|---:|---:|---:|---:|
| arithmetic | 4 / 64 | 1.000 | 0/216 | 1.000 | 39552/43200 (91.6%) | 0 / 0 / 216 |
| arithmetic | 8 / 32 | 0.743 | 194/216 | 0.852 | 39472/43200 (91.4%) | 90 / 109 / 17 |
| arithmetic | 16 / 16 | 0.738 | 171/216 | 0.839 | 37372/43200 (86.5%) | 48 / 160 / 8 |
| geometric | 4 / 64 | 1.000 | 0/216 | 1.000 | 35833/43200 (82.9%) | 0 / 0 / 216 |
| geometric | 8 / 32 | 0.741 | 197/216 | 0.819 | 37346/43200 (86.4%) | 91 / 84 / 41 |
| geometric | 16 / 16 | 0.734 | 174/216 | 0.822 | 32220/43200 (74.6%) | 45 / 160 / 11 |

Available-ratio RMSE is conditional on finite positive point support; it is
not unconditional risk and must be read with all availability/usability
counts in all-comparisons.csv. World/frequency conditions overlap and are
not independent biological replications. No acceptance or coverage guarantee
is inferred from pooled descriptive percentages.

### All views: homogeneous world, arithmetic, both noise laws and frequencies

Each entry aggregates four conditions / 800 simulated calls; model-specific
details and geometric results remain in all-comparisons.csv. This table is
a retrospective descriptive display; complete worlds remain in the CSV.

| View | K8 median variance factor | K4 usable coverage | K8 usable coverage | K16 usable coverage |
|---|---:|---:|---:|---:|
| npc | 0.742 | 757/800 | 754/800 | 561/800 |
| edf10rep1 | 0.555 | 736/800 | 764/800 | 756/800 |
| edf10rep2 | 0.676 | 765/800 | 756/800 | 738/800 |
| fig4untreated | 0.777 | 765/800 | 744/800 | 681/800 |
| rtx_rep1 | 0.729 | 729/800 | 745/800 | 679/800 |
| rtx_rep2 | 0.799 | 759/800 | 767/800 | 748/800 |
| rtx_rep3 | 0.625 | 755/800 | 754/800 | 713/800 |
| sdc_corrected | 0.860 | 629/800 | 652/800 | 606/800 |
| sdc_raw | 0.863 | 635/800 | 609/800 | 586/800 |

## What the result does and does not repair

1. The joint-covariance calculation is now about the actual implemented
   finite-phase moments, including unequal counts and dependence across
   frequencies/windows/angles. It is not a denominator-noise argument.
2. The restricted optimality proof now uses the correct five code-label
   equality orbits. Arbitrary-contribution unbiasedness forces the three
   repeated-position coefficients to agree. It is not a lower bound for
   nonlinear estimators; the N=3 toy counterexample is retained.
3. The generic mathematics is established prior art. The contraction is
   exactly K times the classical third-cumulant polynomial across bins.
   Dependent bins require the source-assignment proof, but do not make
   the polynomial new. PRIOR_WORK.md limits the possible novelty claim.
4. The first-pass population scores combine a sixteen-group code-variance
   estimate with the full retained population’s physical floor. They are
   preliminary descriptive estimates, not exact population optima. The
   finite-budget verification deliberately uses a single consistent law.
5. In particular the EDF10-rep2 prefix can miss rare dominant contributions.
   Its new conditional performance cannot fix the original acquisition’s
   poor Fieller calibration or erase the all-field failure. More codes
   cannot remove a physical/contribution-weight floor.
6. Neither covariance nor hypothetical code allocation validates physical
   code exogeneity, labeling effects, source-specific kernel changes or
   shared errors. Kernel-violation observable targets remain distinct from
   physical attenuation. Randomized MERFISH remains simulation-only.

## Decision

Retain this as a bounded theory/design extension with a concrete conditional
benefit and an explicit failure boundary. It improves explanation and
measurement-design analysis; it does not establish a universal better
estimator, a hardware-ready optimum or closure of the physical-validation
gap. The inspected prior work does not justify marketing this as a new
general cumulant or sketching principle. No manuscript claim upgrade or
replacement of Fieller uncertainty is warranted automatically.

Full derivation: ../THEORY.md. Prior-work comparison: ../PRIOR_WORK.md.
Controls and independent recount: ../test_covariance.py and ../verify.py.
Historical artifacts were read-only throughout.
