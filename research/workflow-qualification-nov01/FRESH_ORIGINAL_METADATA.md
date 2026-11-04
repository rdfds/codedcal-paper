# Original tabular checksum diagnosis, before pixel qualification

The default access API serves Dataverse's ingested TAB representation. Those
complete bodies match the returned lengths and have retained SHA256, but not the
catalogued MD5 for the original CSV. Preserve both facts. Retrieve only the three
original-format files (12555510, 12555511, 12555507) with `?format=original`, each
<=512 KiB, <=1.5 MiB total, 40-second timeout, ten-second global pacing, no retries,
stop on HTTP429/Retry-After and >=1 GiB disk reserve. Compare original body MD5 to
the publisher catalogue before using it. No images or outcomes are requested.
