# Bounded workflow qualification study

Start with [RESULTS.md](RESULTS.md). The completed development matrix does not
demonstrate a readout-aware advantage over the strong mixed conventional suite.
The entire original matrix and all source/compute exclusions remain available.
The separate native source check does not authorize independent comparative
testing after that failed utility gate. No manuscript claim is upgraded.

The reusable command and input/adapter contract are documented in
[workflow_challenges.md](../../tool/docs/workflow_challenges.md). It is separate
from the unchanged estimator command. Native policy selection, reference-first
eligibility, per-policy native gates and readable reports are software features;
they are not evidence of improved qualification accuracy.

Evidence entry points:

- `development-freeze.json`, `development/` and `development-worker.log`: exact
  frozen plan, all native/trial records, complete state and original log.
- `frozen-v1/`: exact implementation pinned before the study, including the older
  API. The worker was not restarted after the shipped API repair.
- `adjudicated-summary.json`, `workflow-decisions.csv`, `recovery-profiles.csv`:
  reference-first replay, all 32 workflow decisions and all 332 challenge pairs.
- `qualification-development.png` and `.svg`: descriptive study figure with
  exclusions and limited replication explicit; not a manuscript claim upgrade.
- `fresh-full-workflow-qualification.json`, `full-workflow-control-calls/`:
  all four controls from the single fixed native hybrid, including the fractional
  catalog's weaker blank-proxy precision. Preserve the raw native-gate flag.
- `portable-full-workflow-input.json`, `fresh-qualified-inputs/input-512.npz`:
  exact processed input retained for reproduction. Original publisher hashes,
  conditional-header failure, exact-range receipts and eight native contexts are
  documented separately. Full TIFF hashes were verified at collection and cannot
  be recomputed from selected retained cells alone.
- `fresh-public-patches/ATTRIBUTION.md`: public source attribution and license.
- `confirmation-decision.json`, `completion-manifest.json` and
  `completion-verification.json`: final decision, file pins and offline closure.

No-fit/no-network checks, from the repository root:

```sh
PYTHONPATH=tool/src:research/workflow-qualification-nov01 python -m pytest -q tool/tests research/workflow-qualification-nov01/test_adjudication.py
python research/workflow-qualification-nov01/verify_completion.py
```

The latter also checks preceding files against the retained baseline Git blobs
and compares pre/post software-repair summaries; use the complete repository
history. Original imported-science verification is `python reproduce/verify.py`.
`analyze_completed.py` replays retained catalogs without fitting; it writes derived
reports, so use a working copy if preserving the exact completion file seal.
The archived `run_development.py` deliberately refuses overwrite/restart and
asserts its original freeze. Do not launch it against the repaired live API.
Fresh execution would require a separately authorized plan, not this completion.
