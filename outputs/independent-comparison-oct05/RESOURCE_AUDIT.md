# Public resources and exact evidence boundaries

The BarDensr author release supplies annotated BARseq crops and an actual
nonnegative-regression decoder. This addresses two weaknesses of our earlier
MERFISH-only comparison: a decoder from a distinct inference family can be
tested in its intended assay, and annotations are not outputs of our reference
decoder. It does not supply physical molecule truth, randomized interventions,
or documented animal independence between its crops.

Source article: [Chen et al., 2021](https://doi.org/10.1371/journal.pcbi.1008256).
The authors already studied demixing under reduced resolution. Repeating that
advantage would not be a new contribution. Our frozen question concerns readout
alignment, separate spatial/label endpoints, and an observable preflight against
both a conventional template workflow and the author model.

The pinned public HDF5 inputs are already local; their full SHA256 hashes and
schema are in `bardensr-schema.json`. Only their schema was inspected before
the initial protocol commit `999b523`; threshold selection uses development
annotations and is explicitly supervised. The hand-curation statuses are
retained. Every predicted interior call counts, including calls not matched to
a good annotation; such calls are not declared molecular false positives.

The separate mouse-colon Dryad release remains unavailable to anonymous download
in the tested interfaces: dataset and file metadata returned 200, download API
401, public file-stream endpoint 403. The archive is 79,277,852,292 bytes, exceeding
available disk and unnecessary for this selective local study. The cultured-cell
Dataverse API returned 403. Access failures are recorded in the corresponding
JSON files; they are not evidence that the datasets or results are invalid. No
authentication bypass, new account, collaborator outreach or paid run was used.

The four reserved TMS kidney fields are not opened here. The prior public
prediction failure and unsuccessful MERFISH regression qualification remain
reported. The distinct BARseq workflow cannot retroactively turn either into
a success or validate the temporal intervention.

Dense evidence tensors are retained in local ignored `.npz` files. Source hashes,
every arm's numerical results, portable primary calls, source settings, failures
and the full tensor/call manifest are committed after verification. Public source
files plus the pinned kernels and scripts regenerate the dense tensors; they are
not discarded to hide unfavorable outcomes. Partial unbounded-template results
are retained separately and excluded from the revised validation analysis.
