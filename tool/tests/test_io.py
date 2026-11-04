"""Unit tests for codedcal.io (Insight3 reader, dax band geometry)."""
import numpy as np
import pytest

from codedcal import io as cio


def _write_insight3(path, rec):
    header = b'M425' + bytes(8) + np.array([len(rec)], '<i4').tobytes()
    path.write_bytes(header + rec.tobytes())


def test_insight3_roundtrip(tmp_path):
    n = 7
    rec = np.zeros(n, cio.INSIGHT3_DTYPE)
    rec['x'] = np.arange(n) + 0.25
    rec['y'] = 100.0 - np.arange(n) * 0.5
    rec['h'] = 1000.0
    rec['frame'] = np.arange(n)
    p = tmp_path / 'toy.bin'
    _write_insight3(p, rec)
    a = cio.read_insight3(p)
    assert a.shape == (n,)
    assert np.allclose(a['x'], rec['x'])
    assert np.array_equal(a['frame'], rec['frame'])
    xy = cio.insight3_xy(p)
    assert xy.shape == (n, 2)
    assert np.allclose(xy[:, 0], rec['x'])
    assert np.allclose(xy[:, 1], rec['y'])


def test_insight3_rejects_bad_magic(tmp_path):
    p = tmp_path / 'bad.bin'
    p.write_bytes(b'XXXX' + bytes(80))
    with pytest.raises(ValueError):
        cio.read_insight3(p)


def test_insight3_rejects_truncated(tmp_path):
    rec = np.zeros(5, cio.INSIGHT3_DTYPE)
    p = tmp_path / 'trunc.bin'
    header = b'M425' + bytes(8) + np.array([50], '<i4').tobytes()
    p.write_bytes(header + rec.tobytes())
    with pytest.raises(ValueError):
        cio.read_insight3(p)


def test_dax_band_range_matches_frame_geometry():
    # tms_fetch.py geometry: 2048x2048 uint16 frames (8388608 bytes),
    # 4096 bytes per row, central band rows 512..1536.
    start, end = cio.dax_band_range(0)
    assert (start, end) == (512 * 4096, 1536 * 4096)
    start, end = cio.dax_band_range(4)
    assert start == 4 * 8388608 + 512 * 4096
    assert end - start == 1024 * 4096
