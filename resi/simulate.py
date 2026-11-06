"""Forward model for stochastically multiplexed localization data (RESI-type)."""
import numpy as np


def complexes(rng, side, density, fractions, spacing):
    """Complex centres under complete spatial randomness; returns subunit positions
    and complex ids. `fractions[m-1]` is the fraction of *complexes* of size m;
    subunits of an m-mer lie on a line (chain) with the given spacing."""
    fractions = np.asarray(fractions, float) / np.sum(fractions)
    sizes = np.arange(1, len(fractions) + 1)
    mean_size = (fractions * sizes).sum()
    n_cx = rng.poisson(density * (side / 1000.0) ** 2 / mean_size)
    m = rng.choice(sizes, n_cx, p=fractions)
    centre = rng.uniform(0, side, (n_cx, 2))
    theta = rng.uniform(0, np.pi, n_cx)
    pos, cid = [], []
    for k in range(n_cx):
        offs = (np.arange(m[k]) - (m[k] - 1) / 2) * spacing
        u = np.array([np.cos(theta[k]), np.sin(theta[k])])
        pos.append(centre[k] + offs[:, None] * u)
        cid.append(np.full(m[k], k))
    if not pos:
        return np.zeros((0, 2)), np.zeros(0, int)
    return np.concatenate(pos), np.concatenate(cid)


def acquire(rng, subunits, complex_id, K, efficiency=0.5, pi=None, link_sd=0.0,
            locs_per_molecule=20.0, round_gain=None, loc_sd=5.0):
    """Label subunits (efficiency, random permanent label), then localize.

    Returns localization positions, labels, molecule ids (labelled subunits) and
    complex ids. `round_gain[c]` scales the expected localizations of label c,
    modelling round-specific imaging (label-dependent counts).
    """
    pi = np.full(K, 1.0 / K) if pi is None else np.asarray(pi, float) / np.sum(pi)
    gain = np.ones(K) if round_gain is None else np.asarray(round_gain, float)
    keep = rng.random(len(subunits)) < efficiency
    site = subunits[keep] + rng.normal(0, link_sd, (keep.sum(), 2))
    cid = complex_id[keep]
    code = rng.choice(K, len(site), p=pi)
    n = rng.poisson(locs_per_molecule * gain[code])
    mol = np.repeat(np.arange(len(site)), n)
    xy = site[mol] + rng.normal(0, loc_sd, (len(mol), 2))
    return dict(xy=xy, label=code[mol], mol=mol, complex=cid[mol],
                site=site, site_label=code, site_complex=cid, counts=n)


def relabel(rng, data, K, pi=None):
    """Redraw molecule labels for a fixed configuration (positions and counts)."""
    pi = np.full(K, 1.0 / K) if pi is None else np.asarray(pi, float) / np.sum(pi)
    code = rng.choice(K, len(data['site']), p=pi)
    return code[data['mol']]
