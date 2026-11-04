# Third-party data and software notices

Project-authored code is distributed under the MIT license in `LICENSE`.
This does not relicense third-party records. Source-derived observations,
codebooks, metadata and retained numerical reductions keep the applicable
source terms. Original identifiers, hashes and transformations are recorded
in the adjacent study manifests and `provenance/import.json`.

| Source | Attribution and original record | Terms |
|---|---|---|
| RESI localization data | Ralf Jungmann; data accompanying Reinhardt et al., Nature 617, 711–716 (2023), https://doi.org/10.1038/s41586-023-05925-9; https://doi.org/10.5281/zenodo.7795826 | CC BY 4.0 |
| SDC-OPR localization data | Cecilia Zaza et al.; https://doi.org/10.1038/s41467-025-60263-w; https://doi.org/10.6084/m9.figshare.28741373.v1 | CC BY 4.0 |
| BARseq3 image data | Huihui Qi, Manjari Anant, Dylan Faltine-Gonzalez and Justus Kebschull; https://doi.org/10.5281/zenodo.19024678; associated preprint https://doi.org/10.64898/2026.05.13.724900 | CC BY 4.0 for the Zenodo deposit; no repository-wide software license inferred |
| BarDensr benchmark and kernels | Jackson Loper and contributors; pinned source `dfdafe9055be322b286369d21e3349204c1f6fc3`, https://github.com/jacksonloper/bardensr; Chen et al., https://doi.org/10.1371/journal.pcbi.1008256 | Repository MIT notice retained in `licenses/bardensr-MIT.txt` |
| STARfinder example data and kernels | Jingyi Ren et al.; https://doi.org/10.5281/zenodo.11176779; https://github.com/wanglab-broad/starfinder at `d6a1e2e0d88a04676d26ad17bf5d20a4ff895161` | Data CC BY 4.0; repository MIT notice retained verbatim in `licenses/starfinder-MIT.txt` |
| 4T1 cell-culture images | Jenkin Tsui and co-depositors; MERFISH4T1Cell-Culture Raw Unaligned Images, https://doi.org/10.7910/DVN/NE0WOP, version1 | CC BY 4.0 |
| 4T1 tissue images | Jenkin Tsui and co-depositors; MERFISH4T1Tissue Raw Unaligned Images, https://doi.org/10.5281/zenodo.15678603 | CC BY 4.0 |
| Whole-cortex BARseq records | Xiaoyin Chen and contributors; perturbation data https://doi.org/10.17632/5xfzcb4kn8.1 and related article https://doi.org/10.1038/s41586-024-07221-6 | Deposit CC BY 4.0; retained article BioC capture states CC BY 4.0 |
| Allen MERSCOPE VA142 regression codebook | Hongkui Zeng and Allen Institute for Brain Science, Mouse whole-brain transcriptomic cell type atlas—MERSCOPE v1:638850 (2023); Yao et al., https://doi.org/10.1038/s41586-023-06812-z; https://api.brainimagelibrary.org/web/view?bildid=ace-den-fix | CC BY-NC 4.0; exact unmodified codebook fixture only, not MIT |
| DART-FISH metadata only | Kian Kalhor, Chien-ju Chen and Kun Zhang; https://doi.org/10.5281/zenodo.8253772 | CC BY 4.0; metadata is not evidence of a completed experiment |

CC BY 4.0 legal terms: https://creativecommons.org/licenses/by/4.0/legalcode.
This project selects and reduces observations, constructs synthetic perturbations,
and computes new statistics; these changes are described in the study protocols.
Attribution does not imply endorsement by any original author.

Dependencies such as NumPy, SciPy, Matplotlib, BarDensr and Serval are installed
separately and retain their own licenses. Calling a native kernel does not mean
that this repository redistributes or exactly reproduces its complete pipeline.
The upstream STARfinder MIT notice retains its original placeholder copyright
line; it has not been silently completed or replaced.

Four nonessential original website/source captures are omitted from this public
snapshot. `provenance/external-source-omissions.json` retains their exact hashes,
sizes, source links and reasons. Complete comparison matrices, reduced numerical
inputs, failures and all supported retained-data replays remain included.

The Allen codebook at `outputs/code-design-sep16/codebook_VA142.csv` is an
unmodified third-party fixture used by the legacy certificate regression test.
Its CC BY-NC 4.0 terms (https://creativecommons.org/licenses/by-nc/4.0/legalcode)
remain separate from the software MIT license. The primary BIL record names
this exact 638850 dataset and its license. No raw image or new atlas outcome
is added. The codebook and the analytic geometric-scope witness were restored
from preserved project records to make the complete shipped test suite runnable.
See `provenance/test-fixtures.json` for exact identities and scope.
