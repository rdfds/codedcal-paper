# Exact cross-label statistics for stochastic multiplexed labelling: analyses

In RESI-type imaging each molecule carries one of K labels (DNA-PAINT docking
sequences), drawn at random and kept for the whole acquisition. Localizations
with different labels therefore come from different molecules. This folder
holds the analyses behind the paper: simulation checks of the estimators, the
public CD20 and Nup96 data of Reinhardt et al., *Nature* 617, 711–716 (2023),
and Nup96 data of Zaza et al., *Nat Commun* 16, 4991 (2025). The estimators
themselves are in the `crosslabel` package (`../crosslabel`).

## Files

| File | Contents |
|---|---|
| `verify_exact.py` | Redraws labels on fixed configurations; Monte Carlo means and variances against the exact expressions |
| `verify_model.py` | Recovers known dimer fractions on simulated fields with the CD20 pipeline unchanged |
| `verify_triplets.py`, `verify_angles.py` | Connected triplets and their angle profiles against true complex geometry |
| `pairs_vs_triplets.py` | Dimer and trimer models with identical pair statistics, separated only by triplets |
| `cd20_members.json`, `fetch_cd20.py` | Byte ranges of the CD20 and NPC members in the public archive; range-request fetcher with CRC-checked extraction |
| `cd20.py` | CD20 loading; package statistics bound to four labels |
| `run_cd20.py`, `run_cd20_molecules.py` | Cross-label correlation of raw round localizations and of published RESI molecule positions |
| `summarize_cd20.py` | Label-pair comparisons, merged-label versions, per-round translations |
| `triplets_cd20.py`, `shape_cd20.py` | Connected triplets and assembly shape (triplet angles) on CD20 molecules |
| `permanence_cd20.py` | Cross-round nearest neighbours at duplicate-molecule distances |
| `registration_cd20.py` | Fiducial-free round registration with a shift-injection check per view |
| `design_cd20.py`, `budget_cd20.py` | Round-length thinning, label merging, label-exchangeability permutation test, round-number design curves |
| `npc.py`, `npc_validation.py` | Nup96 NPCs (32 copies) as ground truth: per-pore factorial moments, corner-level estimates, simulation recovery |
| `remote_zip.py`, `zaza_npc.py` | The same NPC validation on the second laboratory's data, read by range requests |
| `robustness.py` | Sensitivity to radius, baseline, background annulus, block size, NPC grouping and geometry |
| `figures.py` | All figures from saved results |
| `reproduce.sh` | Runs every analysis and figure in order |
| `results/` | JSON summaries and saved pair counts for every view |

## Reproduce

Requires Python ≥ 3.9 with numpy, scipy and h5py, plus `curl` and `unrar`.

```sh
pip install -e ../crosslabel
python fetch_cd20.py ../../resi-cd20.nosync        # CD20 views: ~1.1 GB ranged downloads, ~300 MB kept
python fetch_cd20.py ../../resi-cd20.nosync npc    # Nup96 view
sh reproduce.sh ../../resi-cd20.nosync             # all analyses and figures (about 1.5 h)
```

`zaza_npc.py` reads its 1-MB molecule table directly from the 7.2-GB figshare
archive the first time it runs.

## Main quantities

- **Cross-label g(r):** pair density between localizations (or molecules) of
  different labels, relative to a null in which labels are displaced by 1.5–3 µm.
- **Excess labelled neighbours within 30 nm:** ρ ∫₀³⁰ (g − b) 2πr dr, with b the
  mean of g over 60–100 nm and ρ the pair-weighted molecule density. For a
  monomer/dimer mixture it equals the fraction of subunits in dimers times the
  labelling efficiency, so it cannot exceed the labelling efficiency.
- **Connected triplets (m₂):** mean number of unordered pairs of labelled
  partners within 30 nm per molecule, from three-distinct-label triplets after
  subtracting every term pair correlations produce. Zero for any monomer/dimer
  mixture.
- **Triplet angle profile and linearity index:** angle between a molecule's two
  labelled partners; share of connected triplets above 150°.
- **Assembly stoichiometry:** per-pore cross-label pairs and three-label triples
  give efficiency for known copy number, or copy number and efficiency jointly.
- **Round translations:** per-round shifts from plane-corrected centroids of
  close cross-label displacements, checked by injecting a known shift.
