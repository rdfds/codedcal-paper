"""A single inferred source per coarse group yields deterministic agreement."""
import numpy as np
from sensitivity import code_components, diagonal, phase_sums, Q, ANGLES


def test_singleton_code_diagonal_identity_for_every_code_and_count_pattern():
    rng = np.random.default_rng(26102390)
    for code in range(4):
        for counts in ([1, 2, 3, 4], [7, 3, 11, 5]):
            pan = np.repeat(np.arange(4), counts)
            xy = rng.normal(size=(len(pan), 2))*10+np.array([30, -50])
            z, n = phase_sums(xy, pan)
            zz = np.zeros((1, 4, 4, len(Q), len(ANGLES)), complex)
            nn = np.zeros((1, 4, 4)); zz[0, code] = z; nn[0, code] = n
            assert np.allclose(code_components(zz, nn)[0], diagonal(z, n), atol=1e-10)
