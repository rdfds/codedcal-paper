# Estimator

```sh
python -m pip install './tool[hdf5]'
codedcal smlm-budget events.csv --q 0.10 0.15 --json result.json
codedcal smlm-budget aggregates.npz --uncertainty-variant both --json result.json
```

The CSV header is `group,code,panel,x` with optional positive `weight`.
`code` is the source's permanent physical code, an integer from 0 to 3.
`panel` ranges from 0 to 3: panels 0 and 1 belong to window A; panels 2 and 3
belong to window B. `group` identifies the independent grouping unit.
The coordinate `x` and angular spatial frequency `q` use reciprocal units.
A generic event partition does not implement permanent source coding.

Alternatively, supply an NPZ with `events` columns in that order, or complex
`z` with shape `(groups,4,4,nq)`, real `n` with shape `(groups,4,4)`, and optional
`q` with shape `(nq,)`. Use `codedcal smlm-budget --help` for exact options.

The report retains arithmetic/geometric confidence sets, unbounded or disjoint
support, unsupported outcomes and numerical failures. Exit status and a bounded
legacy variance fit do not establish physical assumptions. The optional
`--legacy-fits` output preserves historical bounded fits for reproduction.
Groups, frequency selection, counts, weights and temporal support must match
the declared target; localizations are not independent acquisitions.

The opt-in `smlm-phase-budget` command retains complex phase and reports its
separate magnitude target and conservative confidence construction. Source-
dependent kernel changes, shared errors, code-dependent selection and physical
motion are not resolved merely by numerical agreement.

The remaining library modules retain the related calibration and digital
challenge implementations used in the comparison records. No universal
calibration, decoder or selector superiority is asserted. See the root methods,
research evidence and reproduction guide for assumptions, tests and limitations.
