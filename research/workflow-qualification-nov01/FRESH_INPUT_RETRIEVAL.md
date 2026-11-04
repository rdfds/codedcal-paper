# First fresh input retrieval plan: header and tabular metadata only

Pin the released Harvard Dataverse DOI 10.7910/DVN/NE0WOP, version 1,
metadata saved in public-metadata/metadata-1.txt. Its file entries are explicitly
unrestricted. This is the cell-culture preparation, not the previously exposed
Zenodo 15678603 tumor release. No comparative result is obtained by this step.

Five sequential requests, no automatic retries, >=10-second pacing, 40-second
request timeout, >=1-GiB free-disk reserve. Each small file body is capped at
512 KiB and its SHA256, publisher MD5 comparison, size, response headers and
representation are retained. Files: 12555510 codebook, 12555511 organization,
12555509 microscope JSON, 12555507 stage positions. Use the standard public
access/datafile endpoint; do not request access-restricted resources or send
messages. A mismatch or unavailable file is preserved before interpretation.

For file 12555727 (merFISH_merged_01_000.tiff), publisher size 232717320 and
MD5 940f522dd601a480acc418de93e308c0, request only bytes 0–65535. Require exact
206, Content-Range, length and retained SHA256. Inspect TIFF metadata only,
without decoding image pixels. This step cannot verify the full TIFF MD5.
Stop network calls on 429/Retry-After and preserve receipts. No full TIFF or
archive is authorized by this small retrieval plan. A later bounded, pinned
image retrieval plan is required before native qualification.
