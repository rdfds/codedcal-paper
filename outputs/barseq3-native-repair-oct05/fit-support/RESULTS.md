# Interior-fit sensitivity

Completed 60/60; failed 0; missing 0.
Full-cohort comparison supported: True.

| Field | Decoder | Original rule F1 | Interior rule F1 | Interior strong round F1 | Interior author model F1 | Harmful selections |
|---|---|---:|---:|---:|---:|---:|
| 0 | correlation | 0.6379 | 0.7941 | 0.5000 | 0.9121 | 0/3 |
| 0 | regression | 0.5649 | 0.6855 | 0.5029 | 0.6892 | 0/3 |
| 1 | correlation | 0.5150 | 0.5150 | 0.5150 | 0.4928 | 0/3 |
| 1 | regression | 0.5040 | 0.5040 | 0.5040 | 0.4500 | 0/3 |

Means use both nonzero cases; zero cases remain in decision/harm diagnostics and full ledgers. Fitting uses only observed interior pixels, with template support checked on every iteration. Smaller fit domains change information as well as padding exposure. Original results remain unchanged. Byte-identical decoder inputs can reuse primary calls, with provenance recorded; they are not independent solver replications. Neither cohort supplies independently known RNA identities or measures temporal T. All three radii, sequence identities, composition and point-proxy/support misses remain per arm.
