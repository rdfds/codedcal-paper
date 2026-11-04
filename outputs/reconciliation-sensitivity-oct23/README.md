# Reconciliation sensitivity: retained evidence and verification

Read INTERPRETATION.md first, then RESULTS.md, all-arms.csv and PROTOCOL.md.
The original 81 settings narrow 68 gaps; five views narrow in every setting.
The post-hoc singleton diagnostic narrows 65/77 supported comparisons and
retains four unsupported. These correlated settings are not independent
acquisitions or physical accuracy measurements.

With NumPy and SciPy installed:

```sh
python outputs/reconciliation-sensitivity-oct23/verify.py
```

This is read-only and uses retained in-repository arrays, not original caches,
Git history, network access or new simulations. `--with-raw` additionally
requires the original localization files at the existing registry locations
and freshly hashes them. Recorded machine-specific input paths are provenance,
not the default verifier's input locations. `--write-receipt` explicitly
refreshes the in-tree receipt; it is for maintainers only.

`groups-*.npz` contains per-group moments, all three quadratic numerator/weight
arrays, source counts and selection coverage. `result-*.json` retains every
grid outcome, baseline replay, complete Fieller topology, conditional paired
interval and leave-one-block-out result. `response-uncertainty.json` contains
all 36 coefficient/observation transformations, keeping old response gates.
`singleton-diagnostic.json` is a separately labeled post-hoc diagnostic.

`dependency-manifest.json`, dependency-revisions.json and
additional-dependencies.json distinguish initial inputs, the transparent
performance-only repair and later supplemental dependency inventory. All
supplemental original files match pre-analysis Git objects. Original packages
and failed tests are unchanged. Tests cover the graph cutoff, boundary-event
membership, signed and unsupported sets, all response fits, rectangle bounds,
cache equality, deterministic singleton identity and read-only verification.

The scientific run used Python 3.10.16, NumPy 1.23.5, SciPy 1.10.1 and h5py
3.10.0; the scoped controls also pass in the separate Python 3.13 audit runtime.
Plotting uses Matplotlib 3.7.5. Run generation scripts only in a disposable
working copy: they are maintainer workflows, not read-only verification.
The sequential runner and per-view logs preserve all completed attempts.

No new D076 images, paid jobs, private data, physical source truth or
randomized-readout intervention are supplied by this package.
