# Reproduce the research

Run commands from the repository root. Use new output locations outside this
checkout, keep the retained evidence unchanged, and run historical research test
suites in separate processes: some use a generic module name `run`.

## Integrity

```sh
python3 reproduce/verify.py
```

This hashes the entire included snapshot, checks the export lineage and rejects
missing/altered files, symlinks and external Git object dependencies. It needs
no packages or network access. The new manifest checks portable copies;
`provenance/import.json` preserves their original source hashes separately.

## Install and test

```sh
python3 -m venv .venv
.venv/bin/python -m pip install './tool[test,hdf5]'
PYTHONPATH=tool/src PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -q -p no:cacheprovider tool/tests/test_smlm.py tool/tests/test_temporal_uncertainty.py tool/tests/test_geometric_uncertainty.py tool/tests/test_weighted_target_contract.py tool/tests/test_temporal_reporting_models.py
PYTHONPATH=tool/src:research/finite-code-covariance PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -q -p no:cacheprovider research/finite-code-covariance/test_covariance.py
```

The selected command above exercises the paper's temporal reporting paths. To
run the entire shipped package regression suite, including historical modes:

```sh
PYTHONPATH=tool/src PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -q -p no:cacheprovider tool/tests
```

Optional native dependencies can cause explicit skips; these are not successful
native-decoder qualification. The two immutable legacy fixtures and their terms
are recorded in `provenance/test-fixtures.json`. No study driver, download or
fit is needed to run the package tests.

The exact package/runtime versions actually used for this export are recorded
in `checks/verification.json`; the commands do not claim all dependency versions
have been tested. The library contains historical modes for reproduction.
The temporal uncertainty commands are described in `tool/README.md`.

## All reference settings and comparisons

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python reproduce/reconciliation.py
python3 reproduce/comparisons.py
```

The first recounts all 81 settings, the singleton-removal diagnostic and 36
reference-response interval images. It uses the original numerical assertions
and tolerances with this snapshot's manifest. Its original workstation-lineage
prologue is preserved in the pinned source revision, not silently treated as
matching the normalized metadata. The second checks complete comparison
cohorts and the negative gain result; it does not re-run fitting or decoding.

## Figure data and all four figures

Figures 1, 2 and 4 use Python 3.13.2, NumPy 2.5.3, SciPy 1.18.1 and Matplotlib 3.11.2. Figure 3 uses the recorded historical environment: Python 3.10, NumPy 1.23.5, SciPy 1.10.1 and Matplotlib 3.7.5. Install the latter from `environments/figures-tested.txt` and pass its interpreter explicitly for an exact PNG replay. Newer Matplotlib versions can change layout while retaining the same plotted values.

```sh
python analysis/prepare.py --output /ABSOLUTE/NEW/data
python analysis/figures.py --data /ABSOLUTE/NEW/data --output /ABSOLUTE/NEW/figures --selection-python /PATH/TO/HISTORICAL/ENV/bin/python
```

Figure 1 is an explicitly illustrative schematic. Figures 2–4 use retained
results; Figure 3 regenerates the complete 81-setting display. No new data,
optimization or outcome study is performed. Compare all six scientific data
exports and the PNG pixels. PDF/SVG creation metadata can differ.

## Complete covariance and interval recount

Use **Python 3.10, NumPy 1.23.5 and SciPy 1.10.1** from
`environments/temporal-exact.txt`. Do not substitute a newer runtime and widen
an endpoint tolerance around a failure.

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/finite-code-covariance/verify.py --prefix research/finite-code-covariance/results-attempt-01 --finite research/finite-code-covariance/finite-budget-attempt-01 --output /ABSOLUTE/NEW/covariance-recount.json
```

This unchanged verifier checks all 126 worlds (including 18 excluded worlds),
7,327 source hashes, 259,200 sufficient-statistic confidence-set cases and 1,296
raw-group checks against the shipped APIs. Original Monte Carlo discrepancies
remain. These are computational checks of the stipulated model, not physical
validation. A new output path prevents accidental replacement of a prior receipt.

The retained study scripts include original generation methods and checkpoint
records. Some write beside their inputs or expect externally obtained originals;
run those only in a disposable independent copy after reading their protocols.
The supported commands above need only included retained observations.

## Retained reporting details

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python analysis/reporting_details.py --root . --output /ABSOLUTE/NEW/reporting
```

This read-only reaggregation reproduces all nine historical/centered-float64
baselines, all 162 spatial-support/validity evaluations, the complete signed
residual categories and all four variant/method gate summaries. It makes no
new fits or random draws and does not change qualification thresholds.
Compare every output byte with `data/reporting/`. The actual spatial resampling
frame uses 7–12 occupied coarse-group tiles of the 4-by-3 partition, retaining
cells with zero selected sources within that frame. The current mathematical
specifications are in `methods/temporal-model.md`, `methods/interpretation.md`
and `methods/covariance-reporting.md`; their numbered literature IDs resolve
in `methods/references.md`. Historical proof/study files remain unchanged.
