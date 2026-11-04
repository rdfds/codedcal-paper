# Preparation-derived temporal measurement validation

9 October 2026. Settings frozen before new measurement-validation outcomes.
Existing exposed public acquisitions only; no new pixels, reserved preparation,
paid/cloud run or physical randomized-readout intervention. This is conditional
simulation validation, not independent empirical confirmation.

## Cohort and simulated observation law

Use all eight distinct acquisitions / nine raw-corrected views in Source Data
69, preserving the Figure4 RTX alias exclusion. Use every candidate event
cluster before source QC, including singletons and unsupported sources. Pool
candidate clusters from the four historical rounds within each coarse group;
these inferred clusters define synthetic emitters, not established physical
identities. No oracle eligibility or code/reference agreement selects inputs.

Retain source geometry, exact event-panel counts, localization multiplicities,
reported event variances and within-panel residuals. New permanent source codes
are independently uniform over four labels. The historical code estimator sees
only aggregate code/panel phase and count sums, never source labels. Full
known-source diagonal Fourier evaluation is the information-rich comparator.
Clustering/QC is an additional sampled reference-coverage diagnostic.

Binding events generate correlated bursts: all localizations in a burst share
its simulated position, with its original multiplicity. This explicitly tests
perfect within-burst correlation, not a full optical/localization pipeline.
Gaussian event noise uses each inferred source's within-panel covariance;
low-support sources use a preparation-wide covariance fallback. The second
model resamples its pooled, centered within-panel residuals, corrected for
lost degrees of freedom, with the same low-support Gaussian fallback. Both
primary laws are stationary per source and independent of new codes/counts.
No Gaussian approximation to finite phase sums is used.

Finite panels retain their original count conditions. Frame boundaries are
reported as provenance; the imposed persistent-window law is not a fitted
physical lag law or reconstruction of real motion. Report candidate geometry,
noise/fallback, all source and nonzero triple-weight counts.

## Worlds and targets

For each view and noise law, primary worlds are: homogeneous T=0 and T=3 nm^2;
group-common motion at total T=3 (half source-specific, half common within a
coarse group); and source-heterogeneous motion with T_s=3 times a normalized
factor increasing with source noise. The factor is fixed from templates before
outcomes. The heterogeneity target is the exact count/kernel-weighted
attenuation at each q, not an arithmetic variance or Gaussian-fit T.

Three separate violations at source-motion T=3 are retained: add 25 nm^2 of
event-noise variance in window B; add source-dependent variance in window B;
and apply static code-dependent x offsets [-15,-5,5,15] nm. The first two test
kernel stability; the third tests code exogeneity. Never treat quiet D as a
certificate. All parameter choices are declared stress settings, not inferred
physical defects in these acquisitions.

Calculate analytic full-population expected cross/reference moments from the
known generator and all templates, averaging over projection angles. Keep
(i) physical displacement attenuation with stationary baseline kernel weights,
and (ii) the observable diagonal-ratio functional under the actual noise law
separate. Kernel violations may preserve extraction while corrupting physical
interpretation; code-exogeneity violations need not preserve extraction.

## Sampling, uncertainty and gates

Primary: 40 independent physical-field realizations per view/model/world, each
with two independent permanent-code assignments. Violations: 12 physical-field
realizations with two assignments. Code assignments sharing a physical field
are not independent biological or Monte Carlo field replicates. Use the
observed number of coarse groups per acquisition as the primary group count;
report a 512-group standardized sensitivity from the same generated draws.
If the observed count is below 512, generate 512 and retain the observed prefix.
Templates are sampled uniformly with replacement, independently by group.
All seeds derive deterministically from view/model/world/physical-replicate IDs.

At q=0.10 and 0.15 nm^-1 use pointwise 95% Fieller sets from independent group
contributions, retaining bounded, disjoint, all-real and empty sets. Coverage
counts unsupported/unbounded sets honestly; usable measurement coverage is
separate and requires a finite positive-ratio interval and positive observed
reference. No bounded optimizer is an uncertainty certificate. Confidence
sets are not adjusted for multiplicity; no familywise claim is made.

Declared conditional qualification, separately per primary world/frequency:
absolute Monte Carlo mean ratio bias <=0.01; >=90% usable calls; >=90% coverage
of the physical attenuation target; and the physical-field cluster-bootstrap
95% coverage interval must include nominal 95%. These are practical study
criteria, not established field standards. Ratio bias is calculated over valid
positive-reference ratios; its availability and potentially extreme tails are
reported beside unconditional coverage. Also retain RMSE and interval widths.
Do not drop failed preparations, change criteria, select a noise model or add
replicates because an outcome is unfavorable.

D is tested pointwise at 5% using group-level log-ratio influence. Record
invalid diagnostics, false alerts in valid worlds, and violations producing
absolute physical-target error >0.01 with a quiet diagnostic. No deployable
recommendation is claimed. Uncertainty summaries resample whole physical fields
(2,000 draws), retaining the paired code assignments.

For homogeneous T=3, sample clustering/QC on the first 16 groups in the first
8 physical realizations, separately by generated code with historical radius
and QC unchanged. Record rejection, assignment impurity, localization coverage
and stationary baseline kernel-weight mass retained by PURE eligible inferred
sources. Mixed clusters cannot be credited with known-source coverage.
This finite diagnostic sample cannot establish the acquisition's full-clustering
runtime or an independent physically verified source-selection rate.

## Execution and reporting

One scientific worker, resumable finite per-world outputs, no competing jobs.
Retain source hashes, source catalogs and generator settings, group contributions,
analytic truth, every case, failed controls, qualification failure and diagnostic
miss. Check exact code-average extraction, analytic Gaussian/residual kernels,
Fieller set membership, persistent-motion identity, and explicit source-label
exclusion before primary outcomes. Small commits and pushes, no merge/submission.
Any implementation repair is recorded; outcome-driven redesign is a separate
study and cannot replace this one. Update manuscript claims only after complete
recount and only where this conditional evidence supports them.

### Implementation details fixed before outcomes

The noise factor is f_s = (0.5 + v_s / median(v)), divided by its arithmetic
mean over all retained inferred sources, with v_s = trace(cov_s)/2. The
source-dependent window-B stress adds 25*f_s nm^2 per axis. Qualification uses
joint usable-and-covering frequency; unconditional Fieller coverage and
coverage conditional on usability are reported separately. The bootstrap gate
also uses joint usable coverage. No infinite confidence set earns a usable call.

The catalog is the complete **retained candidate population** in Source Data
69; upstream event/group construction already excluded insufficient input
groups. This experiment cannot establish the accuracy or coverage of discarded
raw input or prove the inferred clusters are genuine molecules. The z nuisance
is retained event-relative z from the template for sampled QC, not a new 3D
noise validation. Fragmented pure clusters receive only their captured
count-triple fraction of a source's kernel weight.
