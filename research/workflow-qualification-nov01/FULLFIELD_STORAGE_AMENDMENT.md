# Bounded storage repair after a zero-request guard stop

The first context-reader attempt stopped BEFORE network I/O:0 requests,0 bytes,
when available disk fell below the initial1-GiB operational reserve. Its state
and log are retained. No access/rate/integrity failure occurred and no scientific
outcome was observed. The actual maximum added storage is8*3*5,171,328 bytes plus
small NPZ/header overhead (<120 MiB); responses are held in memory, one round at
a time. Require initially512 MiB+120 MiB free, and512 MiB+16 MiB before each body.
This retains a512-MiB safety reserve after the worst-case decoded round write.
Stop again on a failed guard. All scientific thresholds, source pins,32-request/
128-MiB global caps and exact206/ETag/length/hash checks stay in force. Only this
specific zero-request disk stop may be resumed; no blind retries of transfers.
