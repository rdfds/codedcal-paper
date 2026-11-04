# Native preflight timing correction

The first worker's retained log and JSON show BARseq3 field0 regression DID
complete its native control set:141 nonblank calls, passed, cumulative298.51s.
The348-second interruption was during field1 regression, after field1 correlation
also passed. QUALIFICATION_REVISION.md incorrectly described the in-flight work
as the first regression control set; its frozen wording is corrected here.

The compute exclusion stands: one actual field0 regression decode took93.26s,
its unpadded integer-oracle decode89.60s, and fractional oracle112.04s. This is
verified functionality with excessive cost for the bounded new grid, not a
scientific regression failure or an unsupported field0 native decoder. The
new matrix excludes both BARseq3 image-regression challenge grids solely for
that pre-comparison compute reason. All141-call evidence is retained.
