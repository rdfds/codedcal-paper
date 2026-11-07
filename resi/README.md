# Exact cross-label statistics for stochastic multiplexed labelling

In RESI-type imaging each molecule carries one of K labels (DNA-PAINT docking
sequences), drawn at random and kept for the whole acquisition. Localizations
with different labels therefore come from different molecules. This folder
contains estimators that use that fact without assigning localizations to
molecules, checks of their exact properties, and an analysis of the public CD20
data of Reinhardt et al., *Nature* 617, 711–716 (2023).

## Files

| File | Contents |
|---|---|
| `labels.py` | Pair and star-triplet statistics by label; exact same/between-molecule decomposition and label-draw variance |
| `simulate.py` | Forward model: complexes under complete spatial randomness, labelling efficiency, random permanent labels, localization |
| `verify_exact.py` | Redraws labels on fixed configurations and compares Monte Carlo means and variances with the exact expressions |
| `model.py` | Cross-label pair correlation of a monomer/dimer mixture; weighted fit and bootstrap refits |
| `verify_model.py` | Recovers known dimer fractions on simulated fields with the CD20 pipeline unchanged |
| `cd20_members.json`, `fetch_cd20.py` | Byte ranges of the CD20 members in the public archive and a range-request fetcher with CRC-checked extraction |
| `cd20.py` | Block-wise cross-label pair counts, shifted-label null, excess neighbours, registration offsets |
| `run_cd20.py` | Per-view analysis of raw round localizations |
| `run_cd20_molecules.py` | The same statistics on the published RESI molecule positions |
| `summarize_cd20.py` | Label-pair comparisons, merged-label versions and per-round translations from saved counts |
| `verify_triplets.py`, `triplets_cd20.py` | Connected three-label triplet statistic: simulation check and CD20 molecule positions |
| `permanence_cd20.py` | Cross-round nearest neighbours at duplicate-molecule distances |
| `design.py` | Distinct-label probabilities and fixed-budget information for K rounds |
| `design_cd20.py`, `budget_cd20.py` | Round-length thinning, label merging, label-exchangeability permutation test and design curves on CD20 |
| `npc.py`, `npc_validation.py` | Nup96 nuclear pores (32 copies) as ground truth: per-pore factorial moments from cross-label pairs and triples, corner-level estimates, simulation recovery |
| `figures.py` | Main figures from saved results |
| `results/` | JSON summaries and saved pair counts for every view |

## Reproduce

Requires numpy, scipy and h5py, plus `curl` and `unrar` for the fetch step.

```sh
cd resi
python verify_exact.py
python verify_model.py
python fetch_cd20.py ../../resi-cd20.nosync          # ~1.1 GB of ranged downloads, ~300 MB kept
python run_cd20.py ../../resi-cd20.nosync
python summarize_cd20.py ../../resi-cd20.nosync
python run_cd20_molecules.py ../../resi-cd20.nosync
python verify_triplets.py
python triplets_cd20.py ../../resi-cd20.nosync
python permanence_cd20.py ../../resi-cd20.nosync
python design_cd20.py ../../resi-cd20.nosync
python budget_cd20.py ../../resi-cd20.nosync
python fetch_cd20.py ../../resi-cd20.nosync npc
python npc_validation.py ../../resi-cd20.nosync
python figures.py
```

## Main quantities

- **Cross-label g(r):** pair density between localizations (or molecules) of
  different labels, relative to a null in which labels are displaced by 1.5–3 µm.
- **Excess labelled neighbours within 30 nm:** ρ ∫₀³⁰ (g − b) 2πr dr, with b the
  mean of g over 60–100 nm and ρ the pair-weighted molecule density. For a
  monomer/dimer mixture it equals the fraction of subunits in dimers times the
  labelling efficiency, so it cannot exceed the labelling efficiency.
- **Connected triplets (m2):** mean number of unordered pairs of labelled
  partners within 30 nm per molecule, from three-distinct-label triplets after
  subtracting every term pair correlations produce. Zero for any monomer/dimer
  mixture.
- **Round translations:** least-squares per-round shifts from the mean
  displacement of close cross-label pairs, which is zero in expectation for
  correctly registered rounds.
