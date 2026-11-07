# crosslabel

Exact statistics for stochastically multiplexed localization microscopy.

In RESI and related Exchange-PAINT schemes, each molecule is labelled with one of
K orthogonal docking sequences drawn at random, and keeps that label for the
whole acquisition. Localizations (or molecule positions) with **different
labels therefore come from different molecules**. `crosslabel` turns that fact
into exact, clustering-free measures of molecular organization:

| Quantity | What it measures | Key property |
|---|---|---|
| Cross-label g(r) | Pair correlation between molecules | No repeat-localization (blinking) artefact at any K |
| Excess neighbours | Labelled partners within a radius | Model-free; at most the labelling efficiency for monomer/dimer mixtures |
| Connected triplets m₂ | Labelled partner pairs per molecule | Exactly zero for any monomer/dimer mixture; detects higher-order assemblies |
| Triplet angle profile | Angle between a molecule's two partners | Distinguishes chains, rings and compact clusters |
| Label offsets | Registration shifts between labels/rounds | No fiducials needed |
| Assembly stoichiometry | Copy number and efficiency of discrete assemblies (e.g. NPCs) | Needs no resolution of same-label neighbours |
| Design | Distinct-label probabilities, round-number planning | Exact combinatorics for any label frequencies |

## Install

```sh
pip install -e .            # from this directory
pip install -e ".[test]"    # with pytest
```

Requires Python ≥ 3.9, NumPy and SciPy.

## Use

```python
import numpy as np
import crosslabel as cl

# xy: (N, 2) positions in nm (one per molecule, or raw localizations)
# labels: (N,) integers 0..K-1 (docking sequence / imaging round)
rng = np.random.default_rng(0)
res = cl.cross_label_g(xy, labels, rng)                 # g(r) with a shifted-label null
rho = cl.plateau_density(xy)                            # density scale of the normalization
excess, ci = cl.excess_with_ci(res, rho, rng)           # labelled partners within 30 nm
m2, m2_ci, _ = cl.connected_triplets_m2(xy, labels, res, rng)
profile, boot = cl.triplet_angle_profile(xy, labels, res, rng)

# discrete assemblies with known copy number (counts per group and label)
N = cl.group_counts(labels, group_ids)
pi = np.bincount(labels) / len(labels)
cl.assembly_stoichiometry(N, pi, n=32)   # efficiency from counts, pairs, triples; n and p jointly
```

`examples/quickstart.py` simulates two fields with identical pair statistics
(dimers and equilateral trimers) and shows that only the triplet statistics
separate them.

## Interpretation notes

- Labels must be assigned to molecules independently of their position and
  organization, and must not change during the acquisition.
- Excess neighbours and m₂ count *labelled and detected* partners; divide by the
  labelling-and-detection efficiency (or its square) for subunit numbers. A
  co-imaged standard of known stoichiometry (such as Nup96 NPCs) gives that
  efficiency through `assembly_stoichiometry`.
- On raw localizations, background (nonspecific) localizations dilute
  cross-label statistics; molecule positions from per-round clustering avoid
  this. Compare label pairs to check that labels are exchangeable.

## Test

```sh
pytest
```

## License

MIT
