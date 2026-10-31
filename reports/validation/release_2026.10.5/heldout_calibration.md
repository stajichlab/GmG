# Held-out calibration (target precision 0.95)

70791 bench queries (clean labels: 72 suspect genomes dropped), 1909 genus groups, 20 random splits by genus (half the genera calibrate, the other half test).

`held-out precision` = precision of the calls made with thresholds fitted on the other half; brackets = 2.5-97.5 percentile over the splits. `clade-out` = fit on Ascomycota queries, test on all other phyla.


## `known` thresholds

| rank | in-sample precision | held-out precision | held-out coverage | fitted min ANI (mean ± sd) | clade-out precision (coverage, n) |
|---|---|---|---|---|---|
| phylum | 0.965 | 0.962 [0.942, 0.970] | 0.996 [0.994, 0.998] | 0.714 ± 0.006 | 0.921 (0.992, 4165) |
| class | 0.954 | 0.952 [0.920, 0.973] | 0.992 [0.984, 0.996] | 0.714 ± 0.006 | 0.927 (0.990, 3990) |
| order | 0.951 | 0.948 [0.910, 0.978] | 0.985 [0.965, 0.994] | 0.728 ± 0.020 | 0.917 (0.988, 3955) |
| family | 0.950 | 0.944 [0.891, 0.983] | 0.961 [0.883, 0.989] | 0.748 ± 0.033 | 0.920 (0.938, 3770) |
| genus | 0.953 | 0.946 [0.918, 0.972] | 0.980 [0.966, 0.992] | 0.751 ± 0.025 | 0.874 (0.984, 4097) |
| species | 0.950 | 0.949 [0.928, 0.976] | 0.539 [0.321, 0.617] | 0.994 ± 0.002 | 0.893 (0.353, 1502) |

## `novel` thresholds

| rank | in-sample precision | held-out precision | held-out coverage | fitted min ANI (mean ± sd) | clade-out precision (coverage, n) |
|---|---|---|---|---|---|
| phylum | 0.950 | 0.937 [0.864, 0.977] | 0.608 [0.504, 0.709] | 0.791 ± 0.009 | 0.911 (0.446, 1871) |
| class | 0.954 | 0.951 [0.867, 1.000] | 0.204 [0.060, 0.493] | 0.908 ± 0.046 | 0.994 (0.116, 469) |
| order | 0.961 | 0.916 [0.722, 0.996] | 0.138 [0.059, 0.218] | 0.933 ± 0.016 | 0.977 (0.064, 256) |
| family | 0.960 | 0.762 [0.614, 0.896] | 0.045 [0.020, 0.103] | 0.968 ± 0.022 | 0.833 (0.018, 72) |
| genus | 0.950 | 0.937 [0.878, 0.970] | 0.224 [0.172, 0.298] | 0.982 ± 0.009 | 0.891 (0.129, 1033) |
| species | - | 0.819 [0.791, 0.852] | 0.026 [0.020, 0.032] | 1.000 ± 0.000 | - (-, 0) |

Reading: a held-out precision at or above the target means the thresholds generalise to genera they were not fitted on; below it means the in-sample numbers were optimistic. Ranks with no fitted threshold at the target (`-`) never reach it for that pool.

