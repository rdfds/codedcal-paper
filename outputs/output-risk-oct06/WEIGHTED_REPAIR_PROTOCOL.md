# Selection-weighted output sensitivity — post-hoc repair

The first five BARseq3 outcomes of the original profile test are exposed. Its
max-over-compatible statistic gives native false alerts from one-pixel maps
whose mean loss is worse but whose lower block-resampling interval admits them.
Keep that original study unchanged. This separate repair was specified after
those outcomes, before its own weighted forecasts; it is exposed development.

Use the same complete finite map bank, fixed common support, block sufficient
statistics, 500 resamples/seed61021, original incumbent and actual catalogs.
For every resample select the minimum observed loss. Share selection weight
uniformly among maps within1e-6 of that minimum (float32 numerical equality),
so unidentifiable ties do not spuriously select the first array element.
Average incumbent identity disagreement and coding-composition TV using those
selection frequencies. Flag mean identity change >=0.05 or mean TV >0.05.
No outcome labels, known offset, RNA reference, new map class or tuned threshold
enters this forecast. It is expected *selection sensitivity under block
resampling*, not a posterior RNA-error probability or upper coverage bound.

Retain the same min12 blocks/common coverage0.8/incumbent10 coding calls gates.
Compare ordinary/full-bank disagreement, best-fit loss gap, raw fit loss,
decoder confidence and blank/count proxies. Also compare bootstrap-weighted
transform RMS itself: molecular propagation must add value beyond geometrical
uncertainty. All required comparative values are retained without choosing a
winning display threshold. Invisible annotation errors remain untestable.

Numerical controls must catch exact ambiguous-fit false consensus with tie-aware
weights, uniform scaling, missing support and RNA-truth overclaim. No additional
public decoding or unused preparation is authorized by this repair. It cannot
retroactively pass the original gate for the expensive nonnegative extension.
Use the same coverage/sensitivity/false-alert/breadth and strong-comparator
qualification criteria in a separate report. If it fails, retain the failure.
Even a pass requires a separately frozen independent-preparation study.
