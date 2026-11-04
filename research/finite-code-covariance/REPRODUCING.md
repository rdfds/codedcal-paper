# Reproduce the bounded covariance extension

From repository root, after installing NumPy/SciPy and the packaged test extras:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=tool/src:research/finite-code-covariance python -m pytest -q -p no:cacheprovider research/finite-code-covariance/test_covariance.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/finite-code-covariance/run.py --output /ABSOLUTE/NEW/covariance-results
```

The retained source generator was originally run under the environment in
`environments/temporal-exact.txt`; use that environment for exact replay. The
analysis requires all retained Sources 69/70, without original camera data or
private history. It never calls historical mutating drivers. The new output
folder must not already exist, so failed attempts cannot be silently overwritten.

`covariance.py` is a research/design API requiring source-level contributions.
The shipped estimator remains unchanged. It is not a deployable source-blind
standard-error estimator. The 512 code assignments only verify conditional
second moments; they are not independent acquisitions or guaranteed ratio
confidence coverage. The proposed cost functions are stipulated, not measured
hardware costs. Read PROTOCOL.md and THEORY.md before interpreting results.

After the covariance pass completes, reproduce the separately specified
finite-budget verification and independent recount:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/finite-code-covariance/finite_budget.py --prefix /ABSOLUTE/covariance-results --output /ABSOLUTE/NEW/finite-budget
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/finite-code-covariance/verify.py --prefix /ABSOLUTE/covariance-results --finite /ABSOLUTE/finite-budget --output /ABSOLUTE/NEW/verification.json
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/finite-code-covariance/report.py --prefix /ABSOLUTE/covariance-results --finite /ABSOLUTE/finite-budget --output /ABSOLUTE/NEW/report
```

Retained results are in results-attempt-01/ and finite-budget-attempt-01/;
report-final/ contains the complete tabulation. results-final-replay/ repeats
the covariance pass with the final implementation; replay-equality.json checks
that all 126 numerical world files remain byte-identical. Optional orientation,
phase-gradient and invalid-input guards do not change the default calculations.

The independent verifier imports neither covariance.py nor the outcome drivers.
It checks all 7,327 source hashes, reconstructs Gram aggregation/variance scores,
independently recalculates all 259,200 confidence sets from sufficient statistics
and checks 1,296 full raw-group confidence calculations
against the shipped APIs. Its first-experiment raw-group cross-check does not
mean every one of 259,200 cases has an independently retained full group array;
all cases retain means, mean covariance matrices, confidence sets and counts.
Seeds and the retained observation law reproduce the remaining full arrays.

The completed all-interval recount is verification-all-intervals.json. The
initial verification.json is retained as the earlier aggregation/raw-check pass.
The artifact manifest can be verified offline with `python
research/finite-code-covariance/manifest.py --verify`.

Optional descriptive display, using the recorded test environment's Matplotlib:

```sh
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python research/finite-code-covariance/plot.py --prefix /ABSOLUTE/covariance-results --report /ABSOLUTE/report
```

The figure's pooled frequencies are descriptive across overlapping development
worlds; error bars would not provide independent biological uncertainty.
