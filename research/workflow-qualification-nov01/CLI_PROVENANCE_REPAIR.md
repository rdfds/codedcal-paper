# Preserve the paper's historical command source

The curated import manifest pins the original `codedcal/cli.py` as scientific
Source Data. Adding an opt-in subcommand there would invalidate that closure.
Restore its exact baseline bytes and expose the new suite separately as
`python -m codedcal.workflow_cli`. This changes no calibration, challenge plan,
trial or adjudication. The historical estimator and imported numerical source
remain byte-identical. New API/docs/tests are separate, additive files. The
frozen development run never imports the estimator CLI and is unaffected.
