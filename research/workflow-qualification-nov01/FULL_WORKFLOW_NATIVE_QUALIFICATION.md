# Apply native controls to a complete calibration/decoder workflow

The corrected full-plane simple-NN attempt has271 nonblank/39 blank components
and upper95=4.198; it is unqualified. This rejects that particular workflow,
not the public acquisition. A zero-calibration/native-kernel gate is insufficient
for a raw unaligned release when the user supplies a complete calibrated workflow.
The original16-workflow development freeze and all source exclusions stay intact.

Before any fresh comparative outcomes, one strong complete workflow may receive
native controls: actual three-start BarDensr codebook registration (50 iterations;
zero/round/template starts, physical134 gene codewords, named blanks held out of
geometric fitting), followed by the actual Serval CoSiD gain optimizer and native
NN/connected-component extraction. Use exactly10 gain-fit iterations, native
max_dist=.5167, fitting min_norm0/final min_norm1, min_area5/max_area12 for FITTING
only, max_size100000, seed110126, zero entropy/L2 penalties (constructor defaults),
no final area cutoff. Restore source photon units from the retained initial scales
and initialize relative gains as scales/mean(scales), matching the source example's
normalization. This is an explicitly adapted hybrid of strong existing methods,
not a claim of exact complete published Serval reproduction or a new estimator.

Apply the SAME full workflow to native/no-effect/known-inverse/composed-fractional
inputs. Keep50 native calls,99% output support, reference>=10 and blank upper95<=.1.
The single analytically specified512-pixel sampling repair is in NATIVE_SAMPLE_SUPPORT.md; no further crop, threshold, codeword decoding or stopping-criterion changes. This remains
native qualification, not a held-out comparison or independent confirmation.
Run only after the main development worker stops, with a30-minute native-control
limit and retained failures. No parameter grid or replacement run. If it qualifies
and the completed development comparison actually supports a useful decision,
freeze the complete external procedure before any comparative outcomes. Otherwise
report the limitation honestly rather than force a confirmation experiment.
