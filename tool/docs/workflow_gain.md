# Relative gains need explicit amplitude units

For positive per-readout gains `c`, a cosine objective applied to `x/c` is
unchanged by replacing `c` with `lambda*c`, for any positive common multiplier.
An absolute norm cutoff on `x/c` is not unchanged: the norm is divided by
`lambda`. The cosine fit alone therefore does not specify the intensity units
of that cutoff. A numerical gain floor also does not identify physical units.

The opt-in helper declares a reference before outcomes:

```python
from codedcal.workflow_gain import anchor_relative_gains

anchored, diagnostics = anchor_relative_gains(fitted_gains, reference_gains)
```

For positive finite vectors of the same shape it computes

```
log(a) = mean(log(c) - log(reference))
anchored = exp(log(c) - log(a))
```

Thus `mean(log(anchored/reference)) == 0`, and anchoring `lambda*c` gives the
same anchored vector in exact arithmetic. Relative gain ratios and cosine
directions are preserved. The implementation uses logs to avoid products or
ratios overflowing, and rejects an unrepresentable result rather than clipping.

The caller must select and document `reference_gains` and the downstream norm
threshold together. Their units define an explicit inference contract. Changing
from an unanchored gain vector to this reference scale can change historical
calls at the same numeric cutoff. To reproduce the old calls after multiplying
norms by `a`, the old cutoff would also need multiplication by `a`; that would
reproduce the historical convention, not identify a new physical scale.

This helper does not stabilize the relative gain fit, resolve absent readouts,
certify a reference, estimate physical gain or demonstrate molecular recovery.
It is ordinary gauge fixing, not a new calibration estimator. The original
workflow and all frozen results remain unchanged. An unqualified development
adapter and a post-hoc audit of retained failures are in
`research/workflow-failure-audit-nov03`. No old native qualification transfers to
the changed inference contract.
