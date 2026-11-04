# Restore the native conditioning context, before comparative confirmation

The fixed patch's provisional simple-source decoder has186 nonblank and12 blank
components, proxy1.441/upper95=2.335. It does not qualify. Do not call mere native
counts a successful adapter. All prior failures remain immutable.

Inspection identifies a source-context mismatch: native filters and initial
intensity histograms were computed after cropping. The source pipeline filters
whole planes before cropping and learns histograms from whole images. Restore
this context on the SAME field000, middle z-index7 and central256-pixel scoring
patch; no other field, plane, crop, codeword or threshold is selected. Use native
fiducial phase alignment, source highpass2/RL10(window9,sigma1.2)/lowpass1 and
source histogram90th-percentile initial scaling. First qualify native simple
source decoding under the unchanged50-call,99%-support and blank-upper95<=.1
rules. No comparative policy/challenge outcomes may be examined beforehand.

## Finite range retrieval

Only the first8 original TIFF IDs already checksum-verified by whole-file reading
are eligible (RNA bit rounds1..8). Pin every request before I/O. Per file retrieve
first512 bytes and final16384 bytes to parse all45 native IFDs, then exactly the
uncompressed pixel strip of frame6 and the contiguous strips for frames22/23.
Use exact206 Content-Range/length, require the original publisher-MD5 ETag, retain
SHA256 for every range, and require TIFF uint16/1608-square/uncompressed/single-
strip metadata and the expected frame mapping. No pointer guesses or unpinned
fallback reads. Whole-image MD5 remains collection-time verified; this partial
re-read cannot recompute it. Unselected trailing pixel bytes within metadata-tail
ranges are not interpreted. No unselected plane is decoded.

At most32 requests and128 MiB additional transport; each body<=16 MiB, timeout40s,
>=10s global request spacing, no retries, <=30-minute finite worker, exclusive
lock, >=1 GiB disk reserve. Stop on429/Retry-After, ETag/range mismatch, unexpected
TIFF schema, exhausted budget or access error; preserve partial receipts. Retain
full native selected planes only in the local cache and provenance/conditioned
patches for the research evidence. No paid/cloud jobs or new contacts.
