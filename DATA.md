# Data and rights

The repository includes retained reduced localization observations, simulated
fields, sufficient statistics, response controls and complete related comparison
matrices. Source-member names, source identifiers, checksums and selections are
retained alongside each study and in `provenance/import.json`.

Core empirical records are in `outputs/reference-reconciliation-oct07/` and
`outputs/reconciliation-sensitivity-oct23/`. Temporal validation and variants
are in `outputs/temporal-validation-oct09/` and `outputs/temporal-variants-oct17/`.
Every covariance world, excluded case and finite-budget sufficient statistic is
retained in `research/finite-code-covariance/`.

Four historical NPC localization tables are absent from their old cache paths.
The included NPZ observations support the documented retained-data replays;
recreating those reductions requires the cited original release. Full camera
archives and unavailable original cloud artifacts are not reconstructed by
this snapshot. No absence is relabelled as a completed experiment.

Original source publications include Reinhardt et al., *Ångström-resolution
fluorescence microscopy*, Nature 617, 711–716 (2023),
https://doi.org/10.1038/s41586-023-05925-9, with the public release at
https://doi.org/10.5281/zenodo.7795826. Preparation/source-member mappings and
additional dataset identifiers remain in the retained input manifests and
reference-comparison methods; raw and corrected SDC are two views of one acquisition.
The [historical source catalog](provenance/source-catalog.md) and [exact included coverage](provenance/source-coverage.json) identify provenance without claiming this snapshot mirrors every older source package.

Project code is MIT licensed. Third-party data and upstream code retain their
original terms. The source-specific attributions, licenses and transformation
notices are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Published deposit
licenses were checked before release. Four nonessential third-party source/site
captures are omitted, with their identities, sizes and hashes retained in
`provenance/external-source-omissions.json`; the originals remain preserved
separately. No numerical comparison, failure record or supported replay input
is removed. Original camera examples with unresolved redistribution terms are
not included. Access to the preserved development repositories is not required
for the supported retained-data commands in `REPRODUCING.md`.

The randomized molecular-readout intervention remains simulation only. No image experiment from the inaccessible proposed confirmation source or
matched independent physical-motion validation is included.
