# Kernel-robust coherence research extension

Read [REPORT.md](REPORT.md) first. This is an isolated research implementation,
not a physical validation certificate or a replacement for the shipped estimator.

* [PROTOCOL.md](PROTOCOL.md): scope, fixed criteria and allowed calibration data.
* [THEORY.md](THEORY.md): exact bound, covariance, inference and counterexamples.
* [SIMULATION_PLAN.md](SIMULATION_PLAN.md): fixed seed and synthetic controls.
* [PRIOR_WORK.md](PRIOR_WORK.md): narrow novelty boundary.
* `core.py`, `inference.py`: source-coded moments and uncertainty constructions.
* `calibration_feasibility.py`: exposed calibration-only grouping design.
* `simulation_controls.py`: finite power/assumption controls, all outcomes retained.
* `verify.py`, `test_coherence.py`: independent literal code replay and tests.

From repository root, use Python 3.10 with NumPy 1.23.5, SciPy 1.10.1,
NetworkX 2.8.8 for the exact retained calibration optimizer. Tests also passed
in the existing separate Python 3.13 validation environment. Output paths must
be new; the scripts refuse to overwrite retained evidence.

```sh
OPENBLAS_NUM_THREADS=1 python research/kernel-robust-coherence/calibration_feasibility.py --output /tmp/coherence-calibration-replay
OPENBLAS_NUM_THREADS=1 python research/kernel-robust-coherence/simulation_controls.py --output /tmp/coherence-simulation-replay.json
OPENBLAS_NUM_THREADS=1 python -m pytest -q -p no:cacheprovider research/kernel-robust-coherence/test_coherence.py
OPENBLAS_NUM_THREADS=1 python research/kernel-robust-coherence/verify.py --output /tmp/coherence-verification-replay.json
```

The calibration script hashes the existing NPZ and decodes only its five
approved calibration/design arrays. The verifier also reads the frozen
calibration-design JSON to check q=.008, and checks historical artifact hashes;
hashing does not decode scored arrays. The initial verifier expected every
historical manifest to use the same schema and failed on a list-format manifest;
`verification.log` preserves that reader failure. The repaired verifier accepts
both historical schemas without modifying them.

Retained initial and final test logs, exact calibration calculations, all
12,400 synthetic repetitions' aggregate counts/intervals, source/group plan,
and integrity replay accompany the code. The new artifact manifest covers the
stable extension files; it does not claim publisher checksum verification of
external archives or physical independence of any observation.
