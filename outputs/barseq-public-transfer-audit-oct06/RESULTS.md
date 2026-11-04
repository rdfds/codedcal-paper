# Public preparation audit — 6 October 2026

The repaired calibration has a real exposed-data success, but independent
confirmation has not run. The [whole-cortex BARseq paper](https://doi.org/10.1038/s41586-024-07221-6)
lists nine raw acquisitions. Public BIL metadata is accessible; the separate
download host currently times out over HTTPS, HTTP and explicit IPv4 HTTPS.
No image payload was acquired in this audit. Existing files were not deleted.

The release is already known to this project. Source Data 23 / barseq-crosstech-sep21
used D078-1L, cycles 1–4, eight exposed/registered positions, for a spatial-scale
null test and sensitivity controls. Its raw tiles were deleted by its original
runner. That acquisition must not be presented as a fresh holdout. Its result
is spatial cross-technology evidence, not another validated temporal-estimator
modality class. The new raw-image candidates below are different preparations.

| BIL ID | Dataset title/path identity | Specimen identity | Eligibility |
|---|---|---|---|
| ace-dim-pad | D078-1L | D078-2L | Excluded: earlier raw exposure and metadata inconsistency |
| ace-dim-own | D078-2L | D078-1L | Excluded: metadata inconsistency |
| ace-dim-owl | D076-1L | D076-1L | Metadata-qualified new candidate |
| ace-dim-out | D076-4L | D076-4L | Metadata-qualified new candidate |

Candidates are the first two title/directory/specimen-consistent entries in the
paper's listed order after excluding prior raw-image exposure. These choices
use metadata, not molecular outcomes. Their names were searched against this
repository; no earlier candidate raw-image analysis was found. Both are female
C57BL/6J mice, 28 days old, with distinct sample IDs and acquisition dates.
Metadata states four channels ordered G,T,A,C, unlike the BARseq3 adapter's
G,T,C,A; an explicit qualified channel/codebook adapter is necessary.

The dataset paths are:

- https://download.brainimagelibrary.org/23/3e/233eb9bd2df10f2e/20230402_VisReorg_D076_1L/
- https://download.brainimagelibrary.org/23/3e/233eb9bd2df10f2e/20230330_VisReorg_D076_4L/

Metadata also says 104-gene panel, twenty-one optical slices, maximum-projection
gene-sequencing files and an author processing protocol. Those descriptions do
not qualify a decoder adapter or establish RNA truth. Source book, full required
cycles, channel mapping, projection and preprocessing must be audited and frozen
before candidate intensity exposure/scoring. Do not truncate to four cycles or
invent a blank class to force a pass. A metadata-qualified candidate is not a
completed experiment.

The [documented public HTTPS route](https://www.brainimagelibrary.org/download.html)
uses the download host. Globus requires an identity mapping and a support email;
neither was initiated. No account, outreach, private access or paid run was used.
This is an access blocker, not evidence that calibration fails or that a lab
collaboration is scientifically required. Metadata, exact responses and hashes
are retained; access-probes.json records timeouts. A final curl IPv4 HEAD probe
also exited 28 after an 8,006-ms connection timeout and acquired zero bytes.

`TRANSFER_CORE_FREEZE.md` fixes the successful calibration policy and next
sampling/evaluation decisions. Its adapter is explicitly pending, so it is not
misrepresented as a fully frozen independent study or as completed transfer.

## Source-adapter progress after the initial access audit

The author image-processing and MATLAB helper archives have now been downloaded
and verified against publisher SHA-256. A complete seven-round G,T,A,C adapter
is committed in ../barseq-transfer-adapter-oct06/ with original published weights,
instrument profiles and all 111 source sequences (106 coding, five unused).
Eighteen numerical controls pass. All four denoisers load exact published weights;
both actual decoder kernels recover all 106 artificial sources with identity
F1 1 and no blanks. This is synthetic/source qualification, not a real-image
result. Exact original author-runtime/morphology equivalence is not claimed.

Both candidate cycle directories still time out; zero candidate images have
been acquired and none of six planned real native gates have run. The separate
Mendeley data ZIP is processed cell data, not raw sequencing image tensors.
The adapter specification is frozen, while empirical adapter qualification and
the complete independent comparison remain pending. The seven-round source
also invalidates the proposed complete-offset-multiset match between balanced
round and color conditions; the explicit correction in TRANSFER_CORE_FREEZE.md
supersedes that proposal. See ../barseq-transfer-adapter-oct06/RESULTS.md for
executable continuation, declared port differences and retained source hashes.

## Documented service outage identified in the follow-up

The current access pattern matches BIL's announced physical data-center move,
October 5–16, 2026: metadata remains available, while the image filesystem is
unavailable on all systems. The official source was updated October 5, 2026.
October 16 is the operator's expected latest restoration date, not a guarantee;
some services may return earlier. Globus is not an immediate filesystem-outage
workaround. See ../barseq-access-recovery-oct06/RESULTS.md and the
[official notice](https://www.brainimagelibrary.org/datacentermove.html).
Zero candidate raw images and 0/6 real native gates remain. The adapter,
selection rule and numerical-control freeze are unchanged. A finite probe
and optional metadata-resume command are provided in the access-recovery folder.
