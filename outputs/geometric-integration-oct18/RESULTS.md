# Geometric uncertainty integration — 18 October 2026

codedcal 0.4.0 exposes the previously evaluated Source Data 71 geometric
confidence-set construction through the shipped estimator and CLI. It adds no
new physical observations, method-selection rule or empirical validity claim.

## Implemented contract

- `smlm.budget` retains all historical keys and adds `geometric_uncertainty`.
  The existing arithmetic `uncertainty` object and numerical routines are unchanged.
- Each requested level c allocates (1+c)/2 to each window. At c=.95 this is
  two .975 Fieller sets; their positive product/square-root image is the same
  conservative Bonferroni rectangle used in the validated research.
- Window dependence is allowed, independent groups are required, and the marginal
  approximations can fail with concentration or model violations. No simultaneous
  guarantee across frequencies or variants is supplied.
- Complete window confidence sets remain in JSON. Two negative ratios do not
  earn a positive point estimate. A finite transformed subset cannot rescue an
  unbounded, nonpositive or invalid window. Negative Gaussian-equivalent values
  from ratios above one remain visible. This is not uncertainty on a pooled fit.
- `--uncertainty-variant arithmetic|geometric|both` explicitly selects printed
  reports and exit-status checks. Arithmetic remains the default. JSON retains
  both reports and identifies which variants determined status. No automatic
  switch to a favorable estimate occurs.
- The geometric target can differ from arithmetic under heterogeneous motion
  and window weights. Uniform positive-kernel robustness is conditional and can
  sacrifice small-motion resolution; general source/code violations remain.

## Validation

**Exact retained-case replay:** 126 worlds / 3,528 retained fields; all 42,336
geometric cases, all 84,672 window confidence sets and all 42,336 arithmetic
cases match points, complete interval unions, kinds, references and support flags
by exact equality. Every input NPZ and archived case stream is hash-checked.
The Source Data 71 freeze and its unchanged arithmetic Fieller dependency verify.
The five historical fitting functions have identical AST bodies to commit
374e966. No research result file is rewritten by replay.py.

The existing uniform-kernel observed-count coverage counts remain 627/864
arithmetic and 837/864 geometric. Those are retained conditional simulation
outcomes, not new independent observations or improved physical validation.

**196 tests pass.** Analytic controls cover dependent windows, Student rectangle
endpoints, log transforms, positive-set unions/open endpoints, negative and weak
references, degenerate/nonfinite data, confidence allocation and extreme scaling.
Real CLI tests cover CSV/NPZ inputs, JSON, custom levels, unchanged default status,
and supported arithmetic alongside unresolved geometric support. Existing tool,
frozen-replay and completion regressions also pass.

**Built-package smoke passes:** build and install the wheel offline into an
isolated target, then run outside the checkout. Supported CSV with both reports
returns 0; unresolved NPZ geometric returns 2 and retains strict JSON. The wheel
contains version 0.4.0. The first smoke harness resolved the venv executable into
base Python, which lacked NumPy; that harness error and repair are recorded in
wheel-smoke-initial-error.json. No scientific formula changed for that repair.

## Frozen public comparison preservation

The D076 inventory includes all package modules, including temporal modules it
does not use. Its prior full freeze is preserved byte-for-byte as freeze-v2.json,
with source available at 374e966. The pre-image dependency revision updates only
the three changed package module hashes, adds the new geometric module and
preserved freeze, and records the rationale. All spatial policy, endpoint,
adapter and threshold sources remain unchanged. The committed guard and summary
verify: zero images and 0/504 comparison arms. No outcome or exposure justified
this update. The original adapter freeze and Source Data 69/70/71 are unchanged.

## Reproduce

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=tool/src ../validation-venv.nosync/bin/python -m pytest tool/tests tests/test_temporal_frozen_replay.py tests/test_calibration_completion.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 ../merlin-venv.nosync/bin/python outputs/geometric-integration-oct18/replay.py
PYTHONDONTWRITEBYTECODE=1 ../merlin-venv.nosync/bin/python outputs/geometric-integration-oct18/wheel_smoke.py
```

The installed package needs only its declared NumPy/SciPy dependencies. The full
retained replay additionally requires the existing research NPZ/case artifacts;
unit tests for the package itself do not require research downloads. Logs, source
hashes, wheel hash, JSON examples and the freeze revision are retained here.
