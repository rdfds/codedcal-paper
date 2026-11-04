# Native sampling support, derived before complete-workflow qualification

For B=0 named-blank calls, the one-sided95% Poisson upper count is -log(.05).
With134 coding and6 blank classes, the unchanged proxy upper bound is

    U = -log(.05) * (134/6) / N = 66.9046874 / N.

Thus U<=.1 needs AT LEAST670 coding calls even with no blanks. The50-call minimum
is a decoding/recovery reference floor; it does not supply enough proxy precision.
The simple workflow's271 coding/39 blanks also has poor point-proxy quality,
not merely insufficient precision. More pixels alone cannot repair that quality.

Reducing the first patch to256 pixels for compute overlooked this support
requirement. Before any fresh comparative outcome or full calibrated-workflow
native outcome, make ONE sampling correction: same field000, same z7, same
center, same full-plane conditioning, fixed central512 patch. Decoder context
16:496 and scoring/support24:488. No extra field, plane, crop selection, endpoint,
confidence level, minimum count or output threshold is changed. Retain all256
failures. Do not search further sizes if this also fails; report insufficient
support/quality. Run the one complete workflow specified separately, only after
the primary worker stops. This is native design qualification, not comparative
confirmation. The larger input still uses the same four-pair maximum budget per
strategy, with actual runtime/executions/abstentions explicit.
