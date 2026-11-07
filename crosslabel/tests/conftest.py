import numpy as np
import pytest

from crosslabel.simulate import acquire, complexes


@pytest.fixture
def molecule_field():
    """Factory for simulated molecule-level fields (one point per labelled molecule)."""
    def make(fractions, seed=0, shape='chain', side=6000.0, density=150.0, efficiency=0.5, K=4, link_sd=4.0):
        rng = np.random.default_rng(seed)
        sub, cid = complexes(rng, side, density, fractions, 13.5, shape=shape)
        d = acquire(rng, sub, cid, K, efficiency=efficiency, link_sd=link_sd, locs_per_molecule=1, loc_sd=0.0)
        xy = d['site'] + rng.normal(0, 0.5, d['site'].shape)
        return xy, d['site_label'], d['site_complex'], rng
    return make
