# Held-out calibration (target precision 0.95)

> **Re-measured on release 2026.10.5 (2026-10-31, binary `gmg-ca8e3b9`, 72 exclusions): see `release_2026.10.5/README.md` and the files in that directory.** The numbers below are the earlier-release values (BFD bench with 43 suspect genomes dropped); where they differ, `docs/limitations.md` quotes the 2026.10.5 ones.

On release 2026.10.5 the `known` held-out precision is 0.962 / 0.952 / 0.948 / 0.944 / 0.946 / 0.949 (phylum to species) at coverage 0.996 / 0.992 / 0.985 / 0.961 / 0.980 / 0.539: `release_2026.10.5/heldout_calibration.md`.

70878 bench queries (clean labels: 43 suspect genomes dropped), 1914 genus groups, 20 random splits by genus (half the genera calibrate, the other half test).

`held-out precision` = precision of the calls made with thresholds fitted on the other half; brackets = 2.5-97.5 percentile over the splits. `clade-out` = fit on Ascomycota queries, test on all other phyla.


## `known` thresholds

| rank | in-sample precision | held-out precision | held-out coverage | fitted min ANI (mean ± sd) | clade-out precision (coverage, n) |
|---|---|---|---|---|---|
| phylum | 0.965 | 0.962 [0.944, 0.971] | 0.996 [0.994, 0.997] | 0.714 ± 0.006 | 0.921 (0.992, 4171) |
| class | 0.954 | 0.952 [0.922, 0.973] | 0.992 [0.984, 0.996] | 0.716 ± 0.010 | 0.927 (0.990, 3996) |
| order | 0.953 | 0.947 [0.907, 0.978] | 0.982 [0.954, 0.994] | 0.726 ± 0.024 | 0.915 (0.988, 3961) |
| family | 0.950 | 0.945 [0.877, 0.974] | 0.885 [0.398, 0.983] | 0.784 ± 0.057 | 0.945 (0.884, 3559) |
| genus | 0.952 | 0.946 [0.918, 0.972] | 0.978 [0.964, 0.992] | 0.755 ± 0.029 | 0.873 (0.984, 4104) |
| species | 0.950 | 0.948 [0.928, 0.971] | 0.548 [0.466, 0.600] | 0.994 ± 0.002 | 0.892 (0.352, 1500) |

## `novel` thresholds

| rank | in-sample precision | held-out precision | held-out coverage | fitted min ANI (mean ± sd) | clade-out precision (coverage, n) |
|---|---|---|---|---|---|
| phylum | 0.950 | 0.934 [0.864, 0.975] | 0.658 [0.554, 0.745] | 0.787 ± 0.008 | 0.916 (0.456, 1918) |
| class | 0.955 | 0.942 [0.878, 1.000] | 0.279 [0.120, 0.579] | 0.901 ± 0.056 | 0.994 (0.125, 504) |
| order | 0.955 | 0.900 [0.754, 0.996] | 0.150 [0.100, 0.212] | 0.949 ± 0.028 | 0.877 (0.063, 252) |
| family | - | 0.386 [0.345, 0.426] | 0.069 [0.038, 0.119] | 0.973 ± 0.021 | - (-, 0) |
| genus | - | 0.769 [0.651, 0.875] | 0.151 [0.067, 0.213] | 0.993 ± 0.003 | - (-, 0) |
| species | - | 0.810 [0.797, 0.820] | 0.024 [0.020, 0.027] | 1.000 ± 0.000 | - (-, 0) |

Reading: a held-out precision at or above the target means the thresholds generalise to genera they were not fitted on; below it means the in-sample numbers were optimistic. Ranks with no fitted threshold at the target (`-`) never reach it for that pool.

