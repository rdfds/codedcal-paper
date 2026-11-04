# Conditional-request diagnosis before resumption

The32-request context plan first stopped for disk with0 requests, then its first
conditional metadata request returned HTTP412 with0 transported body bytes. The
reader incorrectly used publisher MD5 as an HTTP ETag. The preserved successful
whole-file receipts show multipart ETags (ending `-4`), distinct from whole-file
MD5; the bodies nevertheless matched all publisher MD5s. Pin those actual
collection-time ETags and require exact equality/If-Match on every current range.
This corrects the validator identity, not a source integrity or rate-limit bypass.

Charge the failed request permanently. Reuse the already verified round1 header
prefix (retained65536-byte exact206 response, same collection ETag) rather than
requesting it again. Thus the remaining31 requests plus the failed1 stay within
the unchanged global32-request/128-MiB caps. Preserve all attempt states/logs and
SHA256s. No automatic retry, alternate host, threshold/crop change or new source.
