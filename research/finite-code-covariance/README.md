# Finite-code covariance of temporal calibration

A bounded extension of the permanent-source coded Fourier estimator, developed
27 October 2026. This directory is separate from the preserved manuscript and
numbered historical sources. It does not change the shipped estimator.

- [Results and complete-scope interpretation](report-final/RESULTS.md)
- [Derivation, sampling assumptions and restricted optimality](THEORY.md)
- [Closest prior work and exact classical cumulant relationship](PRIOR_WORK.md)
- [Protocol and dated finite-budget amendment](PROTOCOL.md)
- [Reproduction commands](REPRODUCING.md)
- [Verification record and explicit repairs](CHECKS.md)
- [Independent all-confidence-set recount](verification-all-intervals.json)
- [Final numerical replay equality](replay-equality.json)
- [Complete comparison table](report-final/all-comparisons.csv)
- [Descriptive figure](report-final/covariance-design.png) and [caption](report-final/FIGURE_CAPTION.md)

The conditional covariance decomposes into pair and triple source energies.
After adding the physical-variation floor, it gives an explicit equal-budget
four-versus-eight-code boundary. The finite-state verification finds a median
first-order variance reduction of about 26% at eight codes, with mixed usable
coverage and deterioration at sixteen codes. All model violations and failures
remain visible. The resource model holds counts, photons and kernels fixed by
stipulation; no hardware or physical acquisition improvement is established.

The aggregate polynomial is exactly a scaled classical third-cumulant formula.
The possible contribution is application-specific covariance and measurement
design, not a new general cumulant/sketching principle or universal optimality.
The research API requires source contributions; it is not a source-blind
uncertainty replacement or a deployable workflow recommendation.

ARTIFACTS.json hashes every file in this extension except itself. The source
manifests additionally pin original observations. No original data, source
failure, manuscript claim, packaged numerical module, or original repository
was modified. Run manifest.py --verify for the extension and the repository's
reproduce/verify.py for the preserved scientific snapshot.
