# Final common-scale repair does not rescue these workflows

The precommitted replay completed all 52 decoder predictions. All 26 historical
catalogs reproduced first, and all 24 original challenge assessments matched
exactly before the anchored arm began. Original maps, relative gains, source
scales, native references, support and scoring rules remained fixed.

| Outcome | Barcode model | Group search |
| --- | ---: | ---: |
| Passing integer challenges, historical | 5/12 | 4/12 |
| Passing integer challenges, anchored | 5/12 | 4/12 |
| Passing original held-out challenges, either mode | 1/5 | 1/5 |
| Development workflow rescued | No | No |
| Anchored native calls | 857 | 1,081 |
| Native named-blank proxy upper 95% bound | 0.07807 | 0.06189 |

All 24 identity/composition assessments, native quality measures and both suite
decisions are unchanged. Both suites flag both failing workflows. With zero
adequate examples and one preparation, this remains a tie without an estimable
false-alert or adequate-workflow acceptance comparison. Neither the original
utility criteria nor the failed fractional/native exclusions were relaxed.

The gain vectors actually changed in all 26 cases. A separate post-hoc check of
the retained catalogs found identical class counts and point sets within the
already specified 1e-12-pixel engineering tolerance. The largest matched
coordinate change was 3.44e-13 pixels. Thus the intervention changed the units
but did not change molecular calls in this panel.

The first post-hoc rowwise diagnostic is preserved in posthoc-contrast.json and
contrast.py. It reported only one directly equivalent ordered array, because
near-tied centroid ordinates changed lexicographic row order in 25 cases. Direct
row subtraction is not a valid set comparison: its large coordinate differences
pair different spots. catalog-set-contrast.json instead uses the existing
maximum-cardinality matcher within each class at the unchanged engineering
tolerance; all 26 sets match completely, without unresolved components. These
diagnostics do not alter any original assessment, threshold or primary endpoint.

This rules out the specified final common-scale intervention as a remedy for
these recorded failures. It does not test constrained re-fitting, stabilize
relative gains, or establish the cause of each failure. The earlier retained
audit separately documents representational limits and relative-gain changes;
it remains unchanged. More precise gain units do not validate physical motion.

Close this gain-scale-benefit branch. Keep the explicit amplitude contract as a
software safeguard, without a demonstrated recovery advantage. The focused paper
should retain count matching, supported uncertainty and reference reconciliation
as its substantive results. Do not promote this development test into an
independent confirmation or a workflow-selection benefit.

Resource accounting: 52 completed decoder predictions, four preceding tiny
decoder software controls, one fixed codebook search-index construction, zero
registration/gain optimizations and zero result-cache substitutions. Worker time
was 121.46 seconds (supervisor 128.44 seconds); peak RSS was 423,673,856 bytes.
Per-policy measured decoder unions are in summary.json. They are not separately
measured strategy runtimes or comparable to the cost of full workflow fits.
The conservative reservation is 1,215/1,600; unused old reservations are retained.

No new source calls, physical measurement, fractional study, reserved data,
publication, merge or submission was performed. No further study starts from
this result.
