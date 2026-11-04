# A constructive decoder/codebook limitation

This is a selected algebraic counterexample, frozen before executing the actual
correlation kernel. It is not a prevalence estimate or independent molecular
validation. The public BARseq3 codebook and published correlation settings are
used without changing their threshold (.775) or noise floor (.05).

Two nonnegative, noiseless source words, Kdr (CGACAGTAAA) and Nr2f2
(CTATCTGTGT), each contribute amplitude .5. With known gains, their complete
40-frame signal has a unique nonnegative codebook representation: only those
two code classes fit the exact zero frames, and their columns are independent.
This is an algebraic oracle, not a success claim for a gain-learning solver.

If the G/C and T/A readouts occupy separated loci, the G/C-only signal at one
locus passes the actual correlation threshold for absent Arid1a
(CGCCCGGAGT): evidence 0.77547735. The actual kernel emits this absent class
and no blank class at that point. No true parent passes the threshold. The
complete signal emits no class either: this crowded example is already a
recall limitation for correlation, so it cannot demonstrate successful recovery
by correlation after registration. A nonnegative mixture decoder is a materially
different alternative, not another nearest-vector wrapper.

Thus a zero blank count can coexist with a false gene assignment in this
constructed signal. That does not estimate molecular FDR in a real acquisition
or establish the frequency of such mixtures. The example lies near the fixed
threshold; it is not a robustness claim over arbitrary thresholds/noise/PSFs.
The point-vector test does not include spatial peak selection. All eight blank
scores and every code score were evaluated; no candidate or blank was excluded.

For an isolated signal restricted to k of a ten-readout word's on-frames,
Cauchy–Schwarz bounds its evidence by sqrt(k/10), with the positive noise floor
lowering that bound. At .775, at least seven supported readouts are necessary
even at infinite brightness. In this public codebook, 1,151/1,770 code classes
have no more than six readouts in either camera group. This explains a
conditional loss mechanism, not an assertion that real calls all vanish:
cross-camera PSF tails, overlaps and background violate the restricted-support
assumption. Some blank codes have seven readouts in a camera group, so a
universal claim that all blanks are blocked would be false.

The matched-distribution real-image study is a separate test of whether this
readout structure matters for native-output consistency. Its two patterns have
identical full displacement multisets, covariance and RMS, and it includes
the strong author model and both native decoder families. Only its complete
cohort can support an empirical comparative conclusion. Neither test validates
an unperformed randomized temporal intervention.

Reproduction: run `find_counterexample.py`, retain/freeze the resulting
specification, then run `validate_counterexample.py` with the pinned decoder
environment. Specification, codebook and script hashes are retained in the
JSON records. The actual kernel agrees with the direct algebra to numerical
precision.
