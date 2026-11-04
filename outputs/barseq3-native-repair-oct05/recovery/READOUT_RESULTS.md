# Declared readout-structure diagnostic

Native-reference consistency, not biological sensitivity or true FDR. Groups were fixed from the public codebook before perturbed recovery outcomes.

| Field | Decoder | Error grouping | Policy | Balanced reference recall | Unbalanced reference recall | Blank proxy |
|---|---|---|---|---:|---:|---:|
| 0 | correlation | common | none | 0.0000 | 0.0000 | 0.0000 |
| 0 | correlation | common | common_iterative | 1.0000 | 1.0000 | 0.0000 |
| 0 | correlation | common | population_template | 1.0000 | 1.0000 | 0.0000 |
| 0 | correlation | common | barcode_model | 1.0000 | 0.9677 | 0.0000 |
| 0 | correlation | common | oracle | 1.0000 | 1.0000 | 0.0000 |
| 0 | correlation | color | none | 0.0000 | 0.0000 | 0.0000 |
| 0 | correlation | color | common_iterative | 0.0000 | 0.0000 | 0.0000 |
| 0 | correlation | color | population_template | 0.2308 | 0.0968 | 0.0000 |
| 0 | correlation | color | barcode_model | 0.9231 | 0.9032 | 0.0000 |
| 0 | correlation | color | oracle | 1.0000 | 1.0000 | 0.0000 |
| 0 | regression | common | none | 0.0000 | 0.0000 | 1.5190 |
| 0 | regression | common | common_iterative | 1.0000 | 1.0000 | 0.0000 |
| 0 | regression | common | population_template | 1.0000 | 1.0000 | 0.0000 |
| 0 | regression | common | barcode_model | 0.7805 | 0.8475 | 0.0000 |
| 0 | regression | common | oracle | 1.0000 | 1.0000 | 0.0000 |
| 0 | regression | color | none | 0.0122 | 0.0169 | 0.0000 |
| 0 | regression | color | common_iterative | 0.0122 | 0.0169 | 0.0000 |
| 0 | regression | color | population_template | 0.2073 | 0.0847 | 1.1180 |
| 0 | regression | color | barcode_model | 0.5854 | 0.7288 | 0.0000 |
| 0 | regression | color | oracle | 1.0000 | 1.0000 | 0.0000 |
| 1 | correlation | common | none | 0.0227 | 0.0000 | 0.0000 |
| 1 | correlation | common | common_iterative | 1.0000 | 1.0000 | 0.0000 |
| 1 | correlation | common | population_template | 0.9432 | 0.8378 | 0.0000 |
| 1 | correlation | common | barcode_model | 1.0000 | 1.0000 | 0.0000 |
| 1 | correlation | common | oracle | 1.0000 | 1.0000 | 0.0000 |
| 1 | correlation | color | none | 0.0000 | 0.0541 | 0.0000 |
| 1 | correlation | color | common_iterative | 0.0000 | 0.0541 | 0.0000 |
| 1 | correlation | color | population_template | 0.0114 | 0.0811 | 0.0000 |
| 1 | correlation | color | barcode_model | 0.0114 | 0.0811 | 0.0000 |
| 1 | correlation | color | oracle | 1.0000 | 1.0000 | 0.0000 |
| 1 | regression | common | none | 0.0698 | 0.0000 | 1.9615 |
| 1 | regression | common | common_iterative | 1.0000 | 1.0000 | 0.0000 |
| 1 | regression | common | population_template | 1.0000 | 1.0000 | 0.0000 |
| 1 | regression | common | barcode_model | 1.0000 | 1.0000 | 0.0000 |
| 1 | regression | common | oracle | 1.0000 | 1.0000 | 0.0000 |
| 1 | regression | color | none | 0.0233 | 0.0435 | 0.0000 |
| 1 | regression | color | common_iterative | 0.0233 | 0.0435 | 0.0000 |
| 1 | regression | color | population_template | 0.0233 | 0.1304 | 1.0555 |
| 1 | regression | color | barcode_model | 0.0233 | 0.1304 | 1.0555 |
| 1 | regression | color | oracle | 1.0000 | 1.0000 | 0.0000 |

All 60 arms, including zero cases, all group counts and all three radii are retained in readout-stratification.json. A low or zero blank proxy is not an independent label check. Restricted-support bounds need not describe actual pixels with PSF tails, overlap and background. Synonymous barcodes may cross the camera-balance groups, so gene and exact-sequence recalls are both retained.
