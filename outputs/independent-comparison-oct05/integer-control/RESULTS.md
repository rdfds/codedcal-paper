# Integer-only resampling sensitivity

Post-hoc exposed-crop sensitivity with integer warps and a common prepared decoder window; not independent confirmation.

Completed 120/120 planned arms; failed/missing 0.
Oracle calls exactly reproduce zero/no-correction calls in both families and both crops.

| Crop | Decoder | Strong round-only F1 | Frozen recommendation F1 | Author model F1 | Harmful recommendations |
|---|---|---:|---:|---:|---:|
| development | correlation | 0.7124 | 0.8600 | 0.8640 | 0 / 4 |
| development | regression | 0.6511 | 0.8445 | 0.8440 | 0 / 4 |
| validation | correlation | 0.6363 | 0.7083 | 0.6516 | 1 / 4 |
| validation | regression | 0.6137 | 0.8235 | 0.8326 | 0 / 4 |

Integer application removes interpolation blur; identical oracle calls also control decoder-domain/gain/threshold changes. Applied fitted maps are quantized, so this is a stress test rather than a claim about subpixel performance. The original 312 arms, false recommendation, resampling caveat and external qualification failure remain. All three matching-radius metrics and complete maps/calls are retained per condition.
