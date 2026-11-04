# Existing coded-acquisition reference reconciliation

7 October 2026. Post-hoc validation of exposed acquisitions; this is not a
new independent confirmation or a randomized MERFISH experiment. Historical
results remain unchanged. Settings below precede the new simulation outcomes.

## Cohort and provenance

Include Reinhardt Figure2 NPC, Zaza SDC (raw/corrected views of one acquisition),
EDF10 CD20 replicates 1 and 2, Figure4 untreated CD20, and EDF11/EDF12 RTX CD20
replicates 1, 2 and 3. Figure4 RTX is an alias of EDF11/EDF12 replicate 1;
verify hashes and count it once. No selection by D, frequency agreement or
code/oracle concordance. Eight distinct acquisitions, two source studies;
several preparations belong to the same study/instrument family.

Use historical precision filters, centers, panel boundaries, event linking,
the original clustering radii (10 nm for NPC/CD20, 8 nm for SDC) and source QC. Recover only already-analyzed public
members using the retained range inventory and whole-member SHA-256 identities.
Verify range responses, two extraction implementations and old member hashes.
Keep one raw copy outside Git; remove only newly created redundant copies.
One computationally heavy process at a time. No cloud jobs or paid resources.

## Questions and outputs

1. Reproduce historical oracle source contrasts and counts before interpreting
   a changed reference. Retain mismatches as failures, rather than tuning QC.
2. Measure per-preparation event counts, localization counts, within-panel
   noise, candidate geometry, panel/time support and selection coverage.
   The historical count weight uses *event* counts; code moments use
   *localization* counts. Report both. Oracle means are weighted, code
   Fourier sums are not. Selected and full populations are not interchangeable.
3. Derive the exact four-panel quadratic contrast and distinguish its weighted
   increment product from a Gaussian persistent-window variance and from the
   finite-frequency code attenuation. Do not pool unequal temporal supports as
   a universal lag or regard fitting a common T as establishing its assumptions.
4. Characterize the existing point-A clustering reference with field-derived
   semi-synthetic controls: source geometry and count/noise distributions from
   *all candidate clusters*, before oracle eligibility; known per-window
   Gaussian displacements at T=0, 1, 3, 5 nm^2. Use two independent seeds and
   at least 400 group draws per setting, with the original QC unchanged.
   Report rejected sources, impurity, passing coverage, unweighted/event-count/
   localization-count targets, and reference response. Preserve unsuccessful
   fits. A second residual-noise model tests dependence on Gaussian modeling.
   These controls condition on inferred source assignments and cannot establish
   physical source truth or certify the actual temporal-motion model.
   Also compute an assignment-assisted diagonal Fourier contrast on the exact
   localization rows belonging to eligible clusters. This uses the original
   unweighted Fourier sums, localization counts, frequencies and localization
   time panels, avoiding the event-mean/quadratic-target mismatch. Recompute
   the code contrast on those *same* rows as a conditional diagnostic. This
   oracle-selected restriction is not a deployable assignment-free workflow,
   and selection can break code exogeneity. Retain its original grouped-code population result
   separately; agreement on selected rows does not validate unselected sources.
5. Fit field-specific reference response only when its slope is positive and
   stable across noise models. Before controls run, define the qualification:
   both model slope intervals exclude zero, model slopes differ by at most 10%
   of their mean, and each model's maximum affine residual is <=0.2 nm^2.
   These are declared sensitivity criteria, not established field standards.
   Report model envelopes and raw estimates beside
   any conditional calibrated sensitivity. Never choose a model because it
   moves the reference toward the code estimate. Do not replace a disagreement
   with a fabricated agreement or use interval overlap as an equivalence test.

NPC and SDC existing calibrations are audited alongside the new CD20 controls;
SDC's drift-model mismatch and removed physical-motion intercept must remain
explicit. Repeated controls here characterize the *reference*, not the whole
new estimator's interval coverage or its physical validation.

Execution clarification: the initial EDF10/rep1 400-draw implementation pilot
is retained in pilot400/. It exposes substantial Monte Carlo variation at the
sparse eligible-source rate. Use 4,000 draws uniformly for every primary CD20
control, with the same two seeds, T grid, noise models, QC and qualification
criteria. This increases numerical precision; it does not select a favorable
field, noise model or response threshold. An initial module-name collision
failed before producing any controls; its log is retained and the module now
has an unambiguous name.

Computational clarification: the first three primary fields use the retained
dense clustering implementation. Subsequent controls use sparse radius-graph
connected components with the same strict 10-nm edges, which passes exact
partition controls. This removes quadratic memory allocation without changing
source QC or the simulation law. Paired simulation-group influences retain
common random draws across T within each seed when calculating response-fit
uncertainty; treating those world means as independent would be incorrect.

## Deliverables

Acquisition/alias registry with input hashes; replay and support checks;
per-preparation empirical ingredients; retained per-control outputs; consolidated
comparison table including every failure and disagreement; exact estimand note;
source-ledger/manuscript corrections only where supported. 
