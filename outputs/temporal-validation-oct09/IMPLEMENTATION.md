# Implementation record

Before scientific outcomes, the Fieller control exposed an incorrect set-type
label at exactly zero discriminant with a negative quadratic coefficient.
The confidence inequality then holds for all real ratios, including its
single equality point; it is not a separated complement of an interval.
The type label was repaired. The direct contrast-inequality membership control
and the finite/disjoint/all-real/empty accounting remain mandatory.

Code labels enter only the observation generator. `estimate()` accepts the
coded aggregate phase sums and counts. The comparator alone receives source
phase sums. Noise is generated at event level; burst multiplicity multiplies
phasors and counts, preserving perfect within-event positional correlation.

Retained group arrays have dimensions:
`coded[assignment,group,q,window,(cross,reference)]`,
`known_source[group,q,window,(cross,reference)]`, and sampled template IDs.
All uncertainty and diagnostics can be independently recounted from these
arrays; latent event realizations can be regenerated from the recorded seeds
and hashed templates. Group data are split into small per-field compressed
files; no duplicate event-level raw data are stored.

## Analysis additions after early outcomes (prespecified gates unchanged)

After the NPC/EDF10 replicate-1 fields completed, the descriptive audit showed
that the percentile cluster-bootstrap coverage interval can be degenerate at
all-success coverage. Thus exclusion of nominal 95% from above is not evidence
of undercoverage. Every condition now reports the direction of that failure
and a **post-hoc sensitivity** interval from the bounded-mean Chernoff/KL bound:
`KL(observed field mean || population mean) <= log(2/0.05)/physical_fields`.
Whole-field means lie in [0,1], independently across physical realizations;
no independence between a field's two code assignments is required. This
conservative finite-sample interval never changes the frozen qualification.

To avoid calling a broad finite interval a useful small-motion measurement,
the descriptive analysis also records the fraction of finite positive
intervals that both cover the physical target and exclude ratio=1 (zero
motion). This endpoint was added after early outcomes and is not a frozen
confirmatory test. No generator setting, replicate count, noise choice or
primary gate was changed.

Fieller's contrast inequality is computed exactly for the retained sample
moments, but its nominal Student critical value relies on the joint group-mean
approximation. It is not an exact finite-sample coverage theorem for these
cubic observations. Source and group kernel-weight effective counts and maximum
weight fractions are independently recounted to expose population-weight
concentration. The known-source comparator is also subject to this uncertainty
limitation. A larger number of frozen qualification passes may reflect broader
intervals, not superior measurement resolution.

## Diagnostic-null clarification

Before full-cohort interpretation, the analytic target audit clarified that
source-heterogeneous motion need not give D=1, even under stationary kernels:
window A and B have different source/count/kernel weights. Such alerts are not
false alarms of the common-motion model. Reporting therefore calls only zero,
homogeneous T=3 and marginally source-uniform group-common worlds D-null worlds.
The heterogeneous pooled weighted-attenuation target remains a primary
extraction/uncertainty test. Both baseline and observable source-diagonal
population D targets are retained for every condition. Code-exogeneity stress
can invalidate even the randomized extraction of those diagonal targets.
No simulated outcome, diagnostic value or qualification gate was altered.

The single original scientific worker's execution spans 9–10 October. Settings/generator remained frozen on 9 October;
no replacement worker or scientific rerun was started. Per-field elapsed times
may include inactivity and must not be used as active-compute or method-runtime
comparisons. The execution receipt records actual UTC completion.

The first full recount stopped at its code-provenance check because the
verification script was extended to hash the final execution receipt after
the summary manifest had been created. The failed attempt log is retained.
The manifest is regenerated with final analysis code before recount; no field
observation, truth table, seed or scientific gate changes. This is an analysis
provenance repair, not a failed scientific condition or a replacement run.
