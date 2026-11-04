# codedcal

Count-matched temporal inference and reference reconciliation in localization microscopy.

This research companion contains the estimator, derivations, figure code,
retained observations and complete comparison results supporting
*Count matching clarifies temporal inference and reference agreement in localization microscopy*.

## Start here

- [Install and use the estimator](tool/README.md).
- [Reproduce the results and all four figures](REPRODUCING.md).
- [Current temporal observation model and count-matching proof](methods/temporal-model.md); [retained broader derivations](methods/theory.md).
- [Weighted-input assumptions and bootstrap interpretation](methods/interpretation.md).
- [Conditional coding covariance and design limits](research/finite-code-covariance/THEORY.md).
- [Covariance verification and finite-ratio limits](methods/covariance-reporting.md), [operational diagnostic](methods/diagnostic.md) and [closest prior work](methods/prior-work.md).
- [All baseline, spatial-support and variant-gate details](data/reporting/report.md).
- [Complete reference reconciliation](methods/reference-comparison.md).
- [Data sources and reuse terms](DATA.md), [provenance](PROVENANCE.md) and [verification](checks/README.md).

```sh
python3 reproduce/verify.py
python3 -m venv .venv
.venv/bin/python -m pip install './tool[test,hdf5]'
.venv/bin/codedcal smlm-budget --help
```

The integrity check needs only Python's standard library and reads the included
records. It performs no downloads, fitting or new outcome study. The checkout
includes about 1.4 GB of retained scientific data; allow additional space for
Git objects and the Python environment. Installation of the estimator alone
requires only the `tool/` directory.

## What the method establishes

Under the stated permanent-source coding and observation assumptions, matching
three-panel source-count products removes stable source-specific error-kernel
weighting from a temporal Fourier ratio. The software reports attenuation,
confidence sets, numerical support and conditional variance interpretations.
The covariance analysis retains source-pair/source-triple dependence and an
observation-variation floor.

The empirical comparison retains all eight acquisitions and nine views from
two studies. Matching narrows 68 of 81 correlated selection settings and widens
13; coverage loss and residual disagreement remain. The separate digital
qualification comparisons retain their losses, ties, exclusions and negative
gain replay. Neither conditional inference nor selected-reference agreement
establishes independent physical-motion accuracy.

## Repository layout

| Directory | Contents |
|---|---|
| `tool/` | Installable estimator and tests |
| `analysis/` | Figure-data export and plotting scripts |
| `data/figure-values/` | Complete machine-readable main-figure values |
| `figures/` | Four main displays as PDF, SVG and PNG |
| `methods/` | Observation model, proofs and complete temporal/reference tables |
| `outputs/` | Retained reduced observations, simulations, controls and responses |
| `research/` | Covariance and complete related comparison matrices |
| `reproduce/` | Read-only integrity and numerical recount commands |
| `provenance/` | Original source identifiers, hashes and curation lineage |
| `checks/` | Verification summary and its limits |

The manuscript, cover letter and author declarations are separate submission
documents. Project code is MIT licensed; third-party data retain their own terms.

## Version and source terms

This public repository starts with one initial scientific snapshot. Development
history is maintained separately; no earlier commit ancestry is published here.
Reproduction uses the included files and their content hashes. Original source
software revisions remain identified for scientific provenance.
See [third-party notices](THIRD_PARTY_NOTICES.md) and the exact
[source-capture omissions](provenance/external-source-omissions.json).
