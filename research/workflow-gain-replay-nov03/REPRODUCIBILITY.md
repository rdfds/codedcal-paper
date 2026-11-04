# Reproduce the evidence without another decoder execution

Run from the repository root using the retained production interpreter:

```sh
PYTHONDONTWRITEBYTECODE=1 ../decoder-venv.nosync/bin/python research/workflow-gain-replay-nov03/verify.py
PYTHONDONTWRITEBYTECODE=1 ../decoder-venv.nosync/bin/python research/workflow-gain-replay-nov03/contrast.py
PYTHONDONTWRITEBYTECODE=1 ../decoder-venv.nosync/bin/python research/workflow-gain-replay-nov03/catalog_set_contrast.py
python3 -B reproduce/verify.py
```

The first command verifies pinned inference code, prior seals, all saved maps and
gains, original native references, re-scored assessments, all 26 historical
reproductions and every anchored gain transformation. It replays the frozen
adjudication and checks stable outputs. It does not decode an image or fit a
model. The two contrast commands read retained catalogs only and require the
same results on repetition. The last command checks the 19,089 original imported
scientific files and the numbered-source ledger. The preceding full check is
saved in preceding-reproduction-check.log.

freeze.json records Python 3.10.16; NumPy 1.26.4, SciPy 1.15.3, pandas 2.3.3,
scikit-image 0.25.2 and scikit-learn 1.7.2, plus the existing TensorFlow 2.16.2
package version. The fixed decoder runs no registration or gain optimization;
no TensorFlow calibration is invoked. The entire pinned upstream Serval Python
source is checked.

prepare.py reconstructed every old integer input and verified its completed-call
fingerprint before the freeze. It did not predict outcomes. run.py required that
freeze to be committed and ran through the existing finite process supervisor.
historical-gate.json records the complete first-stage gate. run/ contains every
native/challenge catalog and assessment in both modes. worker.log and
supervisor-state.json preserve the completed finite process. No result was
substituted from a cache.

The four software checks are recorded in software-tests-repaired.log, with four
tiny actual source-decoder predictions in software-control-receipt.json. The
initial import/collection failure is retained in software-tests.log; it made no
decoder predictions. The software-control test and run.py deliberately refuse
to overwrite or silently restart their recorded executions. Do not delete their
receipts to rerun this bounded study. Pure gate checks can be repeated with:

```sh
PYTHONDONTWRITEBYTECODE=1 ../decoder-venv.nosync/bin/python -m pytest -q -p no:cacheprovider research/workflow-gain-replay-nov03/test_replay.py -k 'not actual_source_predictions'
```

INTERPRETATION.md distinguishes the frozen rescue result from the subsequent
ordered-array and set-comparison diagnostics. Both versions are retained. The
engineering comparison tolerance does not replace or loosen molecular scoring.

The 192 files in the three immediately preceding seals are verified unchanged,
including all earlier complete comparisons, post-hoc failures, native controls,
and fractional exclusions. Frozen earlier artifacts were not edited. New work
is additive and exposed development on one public preparation, not an externally
registered design, independent confirmation, molecular truth or physical
validation. The original estimator and all historical numerical modules remain
unchanged.
