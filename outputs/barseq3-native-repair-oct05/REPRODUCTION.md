# Reproduction and interpretation boundaries

Use one scientific process at a time. See PROTOCOL.md, the pinned release
metadata and environment.json. The complete input is a single 934,285,696-byte
public HDF5. Neither AnnData expression values nor private/cell-type labels
are inputs to this study. Source indices are used because the file’s list
index to published FOV-name mapping is unverified.

1. Run `fetch_source.py` and `fetch.py` with the analysis Python. They verify
   the immutable author source SHA256 and the published whole-file MD5.
2. Run `condition.py controls`, then `check_background.py` with the decoder
   Python. The tests compare actual dependency behavior and crop-restricted
   resampling, and exercise the native partial-affine translation sign.
3. Run `condition.py geometry --index 0`, then `condition.py profile`.
   Commit the normalization profile and implementation before new decode
   outcomes. Preserve every per-round geometry, fit diagnostic and source hash.
4. Run `condition.py condition --index 0`; then `qualify.py --index 0`
   separately with `--decoder correlation` and `--decoder regression`.
5. Only when both source-zero gates pass, process fixed source-one geometry,
   conditioning and the two qualifications. Source one’s scoring ROI is the
   previously fixed [1408,1536) spatial holdout. Do not move it or tune cuts.

All scripts run from the repository root using full script paths. The finite
`execute.py` driver performs these stages sequentially, commits/pushes retained
results, stops a failed prerequisite and writes execution-state.json. It also
runs the separately specified INTEGER_CONTROL.md sensitivity. Do not launch
another executor while its recorded PID matches its command. No driver step
reruns a cloud job or launches a replacement scientific design.

Geometry retains the original uint16 assignment casts and published descriptor
retention, RANSAC, spectral inverse, common mask and crop rules. The pinned
BarDensr background implementation applies only its last-axis blur. That
behavior is preserved and tested; the release does not pin its installed
dependency version. One source-zero field supplies normalization percentiles,
where the author configuration uses five unreleased calibration fields.
Only final resampling output is restricted to the ROI; the registration fits
use the full native field. Camera and round fits have no translation fallback.

The correlation kernel and cut are published defaults. The 100-step regression
family is an explicitly added comparator, not the release’s default kernel.
Native gates use point blank-adjusted proxies, not independently known FDR.
Only eight blank classes are present, so even zero observed blanks cannot
certify a low molecular error rate. Shared sequence classes remain pooled.

Passing native gates qualifies an adapter. It does not demonstrate a calibration
benefit, a successful decision rule, independent prospective validation or a
performed randomized-readout intervention. The earlier failed partial adapter
remains in ../barseq3-transfer-oct05. Any subsequent recovery study must be
frozen separately and include the strong author model, failures and abstentions.

## Fixed downstream comparison and verification

The separately committed RECOVERY_PROTOCOL.md and recovery.py define sixty
matched-distribution arms. The finite recovery_execute.py waits for the native
and integer-control prerequisites, commits each field, summarizes the complete
cohort and refreshes the local delivery. Do not change its frozen numerical
script after outcomes. Full-offset equality does not imply equality of temporal
T; read READOUT_STRUCTURE_ANALYSIS.md and PRIOR_ART_BOUNDARIES.md.

FIT_SUPPORT_SENSITIVITY.md defines a post-hoc sixty-arm check of empty shift
borders. Its fit_support_execute.py waits for the primary driver and verifies
completion before running one field at a time. It retains failed fits and
exact input reuse. It does not provide another independent preparation or
replace the primary results.

After both finite drivers finish, run these scripts with the analysis Python:
readout_analysis.py, runtime_analysis.py, verify_results.py,
verify_correlation.py, verify_fit_support.py, verify_assignment.py and
plot_recovery.py. Their reports distinguish retained-call recounts, independent
numerical correlation checks, global matching verification, descriptive runtime
accounting and biological inference. Inspect the rendered comparison figure.
All sixty planned outcomes in each cohort must appear in its ledger, including
failures; no partial cohort supports a complete comparative claim. Solver-input
reuse is verified byte-for-byte and is not another numerical replication.
