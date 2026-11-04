# Uncertainty scope and finite-field sensitivity

Primary uncertainty is the frozen nominal Fieller construction from independent
coarse-group contributions, with a Student critical value and group covariance.
Its polynomial membership test is algebraically exact; the critical-value
calibration relies on a joint group-mean approximation. Population-weight
concentration can undermine that approximation. Qualification is evaluated,
not inferred from a passing formula check.

Coverage is scored against the full template-population physical target. Each
physical-field realization has two new code assignments sharing event noise
and motion. Its mean joint usable-coverage score Y is in [0,1]. Bootstrap draws
resample these independent physical-field scores, preserving code dependence.
At all-success coverage its percentile interval is degenerate. Exclusion of
nominal 95% from above does not demonstrate undercoverage. Frozen gates remain
unchanged, with the failure direction retained.

## Post-hoc bounded-field interval (not a substituted qualification gate)

For independent identically distributed field scores with E[Y]=p and Y in
[0,1], convexity gives E[exp(tY)] <= 1-p+p exp(t). Multiplying these bounds
and minimizing the exponential tail bound gives

P(mean(Y) >= y) <= exp(-N KL(y || p)), for y>p,

and the corresponding lower-tail bound for y<p. Inverting both one-sided
bounds with tail probability alpha/2 yields the confidence set

KL(observed mean || p) <= log(2/alpha)/N.

Here KL(y || p)=y log(y/p)+(1-y) log((1-y)/(1-p)), using its continuous endpoint
limits. A union bound gives at least 1-alpha coverage of the field-score mean,
without independence between the two assignments in a field. At all-success
coverage its lower endpoint is exp(-log(2/alpha)/N), rather than one.

This conservative finite-field interval concerns the repeated-simulation
**coverage probability** under the conditioned generator. It is not a confidence
interval for physical displacement, does not establish biological replication,
and does not remedy a biased estimator. The calculation was added after early
outcomes and is reported for every condition, without altering frozen verdicts.

## Resolution sensitivity (post hoc)

A usable Fieller call requires a finite entirely positive-ratio interval and a
positive observed reference. That condition does not require a sufficiently
narrow interval. The report additionally counts usable intervals that cover
the target and have upper endpoint <1, excluding zero motion in the direction
of positive displacement attenuation. No physical-range clipping is applied.
These pointwise descriptive fractions have no multiplicity-adjusted or frozen
confirmatory interpretation.
