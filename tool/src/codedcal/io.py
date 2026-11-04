"""I/O for coded localization data.

- Insight3 ``.bin`` molecule-list reader (STORM/MERFISH-era format,
  ``M425`` magic, 72-byte records).
- Anonymous HTTP ranged fetch of one band of one frame from a remote
  ``.dax`` movie (raw uint16 little-endian frames), with SHA-256 provenance.
- Optional HDF5 localization-table reader (Picasso-style ``locs`` dataset;
  requires the ``codedcal[hdf5]`` extra).

Provenance (formats and numerics copied verbatim):
  INSIGHT3_DTYPE, read_insight3 <- outputs/codedfish-port-sep16/real_chen2015.py (REC, read_i3)
  fetch_range, dax band geometry <- outputs/codedfish-port-sep16/tms/tms_fetch.py
"""
import hashlib
import time
import urllib.request
from pathlib import Path

import numpy as np

# Insight3 record layout (real_chen2015.py REC, verbatim). 72 bytes/record.
INSIGHT3_DTYPE = np.dtype([('x', '<f4'), ('y', '<f4'), ('xc', '<f4'), ('yc', '<f4'),
                           ('h', '<f4'), ('a', '<f4'), ('w', '<f4'), ('phi', '<f4'),
                           ('ax', '<f4'), ('bg', '<f4'), ('i', '<f4'), ('c', '<i4'),
                           ('density', '<i4'), ('frame', '<i4'), ('length', '<i4'),
                           ('link', '<i4'), ('z', '<f4'), ('zc', '<f4')])


def read_insight3(path):
    """Read an Insight3 .bin molecule list; returns the structured array.

    From real_chen2015.py read_i3 (asserts converted to exceptions; parsing
    identical: 'M425' magic, int32 record count at byte 12, records at 16).
    """
    b = Path(path).read_bytes()
    if b[:4] != b'M425':
        raise ValueError(f'not an Insight3 M425 file: {path}')
    n = int(np.frombuffer(b, '<i4', 1, 12)[0])
    if 16 + 72 * n > len(b):
        raise ValueError(f'truncated Insight3 file: {path} ({n} records, {len(b)} bytes)')
    return np.frombuffer(b, INSIGHT3_DTYPE, n, 16)


def insight3_xy(path):
    """(n, 2) float array of raw x, y from an Insight3 .bin list.

    From real_chen2015.py read_i3 return value (verbatim column stack).
    """
    a = read_insight3(path)
    return np.column_stack([a['x'], a['y']]).astype(float)


# --------------------------------------------------------------------------
# Ranged .dax band fetch (tms_fetch.py; anonymous GETs, raw bytes held in
# memory only, SHA-256 of every chunk recorded).

def fetch_range(url, start, end, retries=4):
    # From tms_fetch.py (verbatim).
    """Ranged GET [start, end) -> bytes. Retries with backoff."""
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end-1}"})
            with urllib.request.urlopen(req, timeout=180) as r:
                data = r.read()
            if len(data) != end - start:
                raise IOError(f"short read {len(data)} != {end-start}")
            return data
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2.0 * (attempt + 1))
    raise IOError(f"failed after {retries} attempts: {url} [{start},{end}): {last}")


def dax_band_range(frame, row0=512, row1=1536, frame_h=2048, frame_w=2048):
    """Byte range [start, end) of rows [row0, row1) of one frame of a raw
    .dax movie (uint16 LE, frame k at bytes [k*frame_bytes, (k+1)*frame_bytes)).

    From tms_fetch.py band_range/frame geometry (verbatim; the constants
    2048x2048 / rows 512..1536 are that experiment's defaults).
    """
    frame_bytes = frame_h * frame_w * 2
    row_bytes = frame_w * 2
    off = frame * frame_bytes
    return off + row0 * row_bytes, off + row1 * row_bytes


def fetch_dax_band(url, frame, row0=512, row1=1536, frame_h=2048, frame_w=2048,
                   retries=4):
    """Fetch one band of one .dax frame; returns (float32 image, manifest).

    The manifest records url, frame, byte range, rows and the SHA-256 of the
    raw chunk, matching the tms_fetch.py get_band provenance record.
    """
    start, end = dax_band_range(frame, row0, row1, frame_h, frame_w)
    raw = fetch_range(url, start, end, retries=retries)
    manifest = {"url": url, "frame": frame, "byte_range": [start, end],
                "rows": [row0, row1], "n_bytes": len(raw),
                "sha256": hashlib.sha256(raw).hexdigest()}
    img = np.frombuffer(raw, dtype="<u2").reshape(row1 - row0, frame_w)
    return img.astype(np.float32), manifest


# --------------------------------------------------------------------------
# Optional HDF5 localization tables (e.g. Picasso/RESI exports).

def read_hdf5_locs(path, dataset='locs'):
    """Read a localization table from an HDF5 file (structured array).

    Requires h5py: pip install 'codedcal[hdf5]'.
    """
    try:
        import h5py
    except ImportError as e:
        raise ImportError("read_hdf5_locs requires h5py; "
                          "install with: pip install 'codedcal[hdf5]'") from e
    with h5py.File(path, 'r') as f:
        if dataset not in f:
            raise KeyError(f'dataset {dataset!r} not in {path}; '
                           f'available: {list(f.keys())}')
        return f[dataset][()]
