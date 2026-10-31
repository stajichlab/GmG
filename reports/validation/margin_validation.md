# Safety margin at genus and family (#58)

Decision record. Date 2026-10-30. The maintainer decided to apply a calibration safety margin at the **genus and family ranks only**; phylum, class, order and
species keep their cutoffs (a margin collapses species coverage 0.55 to 0.11, `clade_calibration.md`). The cutoffs were derived with the held-out protocol of
`clade_calibration.md`, not chosen by hand.

## What changed

| rank | before (ANI, vote consistency) | after | source |
|---|---|---|---|
| family | 0.7715, 0.9 | **0.8837, 0.9** | fit at precision 0.97 on the full `known` pool (ANI 0.8837, consistency 0.5), consistency kept at the shipped 0.9 |
| genus | 0.90, 0.5 | **0.90, 0.9** | the 0.97 fit is ANI 0.8033, consistency 0.9: below the shipped ANI 0.90 (raised in #57), so only the consistency tightens |
| phylum, class, order, species | unchanged | unchanged | |

Rule (`scripts/margin_calibration.py`): per rank the new cutoff is the stricter, component by component, of the shipped cutoff and the fit at 0.95 + 0.02.
A margin may only tighten: the balanced genus ANI of 0.90 was raised on purpose to cut false names for absent genera (`novelty_curve.md`) and is not lowered to a fit of 0.80.
Only the thresholds changed: `Thresholds::default()` (crates/gmg-core/src/call.rs), `data/thresholds/bfd_s10000.known.json` and `.balanced.json`.
`permissive` and `novel` are untouched. **No index rebuild and no release asset is needed**: thresholds are compiled into the binary (and the WASM/Python wrappers), they are
not stored in the GMGIDX03 index. A new binary release carries the change; `tentative_genus_min_ani` (0.7765) is unchanged.

## Derivation (`margin_derivation.md`, `.tsv`; `scripts/hpcc/35_margin_calibration.sbatch`)

70,878 bench queries (43 suspect genomes dropped), 10 random genus splits (half the genera fit, the other half score), pool `known`, target 0.95, margin +0.02.
Held-out precision (called share). `shipped` = the cutoffs before the change; `fit` = refit at 0.95 (the old calibration); `margin` = refit at 0.97; `new` = what ships.

| rank | shipped | fit 0.95 | margin fit 0.97 | **new** |
|---|---|---|---|---|
| family, all | 0.939 (0.94) | 0.929 (0.96) | 0.946 (0.82) | **0.947 (0.82)** |
| genus, all | 0.979 (0.88) | 0.944 (0.98) | 0.965 (0.94) | **0.982 (0.87)** |
| family, Basidiomycota | 0.965 (0.89) | 0.947 (0.93) | 0.974 (0.75) | 0.975 (0.75) |
| genus, Basidiomycota | 0.955 (0.78) | 0.910 (0.96) | 0.937 (0.89) | **0.960 (0.77)** |
| genus, other phyla | 0.955 (0.78) | 0.860 (0.96) | 0.915 (0.88) | 0.958 (0.77) |
| family, other phyla | 0.899 (0.88) | 0.867 (0.92) | 0.905 (0.79) | 0.915 (0.78) |

**The issue's genus numbers (0.944 to 0.965; Basidiomycota 0.910 to 0.937) compare the 0.95 fit with the 0.97 fit. The shipped balanced default (#57) was already past
both** (genus 0.979, Basidiomycota 0.955), so on genus the margin adds +0.3 points (0.979 to 0.982; Basidiomycota +0.5) for one point of coverage. The two points the issue asked for
had already been delivered by #57. Family is where the margin acts: +0.8 points of held-out precision for 12 points of coverage.

## Validation: before and after on real data

Same binary (`gmg-a9acf5c`) and index (`bfd-2026.10.3-s10000.gmgi`) before and after; only `--thresholds` differs (`scripts/hpcc/23_truth_run.sbatch THF=`, `25_release_stability.sbatch THF= TAG=`,
`36_margin_external.sbatch`), scored by `37_margin_score.sbatch` (`validate_run.py`; truth-set calls with decoy calls blanked as in `truth_set.md`) and `scripts/margin_tables.py`.
Phyla other than Ascomycota and Basidiomycota are pooled as "other phyla". Ranks other than family and genus are identical before and after, as intended.


## Post-freeze (1,257 assemblies not in BFD)

Cells: called share / precision [95% Wilson CI] (n with truth). The CI is on the calls made.


**all**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.946 / 0.990 [0.982, 0.994] (n=1257) | 0.946 / 0.990 [0.982, 0.994] (n=1257) | +0.000 | +0.000 |
| class | 0.944 / 0.991 [0.983, 0.995] (n=1238) | 0.944 / 0.991 [0.983, 0.995] (n=1238) | +0.000 | +0.000 |
| order | 0.941 / 0.986 [0.978, 0.992] (n=1238) | 0.941 / 0.986 [0.978, 0.992] (n=1238) | +0.000 | +0.000 |
| family | 0.906 / 0.984 [0.975, 0.990] (n=1250) | 0.850 / 0.992 [0.984, 0.996] (n=1250) | +0.007 | -0.057 |
| genus | 0.841 / 0.988 [0.979, 0.993] (n=1240) | 0.839 / 0.988 [0.979, 0.993] (n=1240) | -0.000 | -0.002 |
| species | 0.750 / 1.000 [0.910, 1.000] (n=52) | 0.750 / 1.000 [0.910, 1.000] (n=52) | +0.000 | +0.000 |

**Ascomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.963 / 0.992 [0.984, 0.996] (n=1041) | 0.963 / 0.992 [0.984, 0.996] (n=1041) | +0.000 | +0.000 |
| class | 0.962 / 0.991 [0.983, 0.995] (n=1041) | 0.962 / 0.991 [0.983, 0.995] (n=1041) | +0.000 | +0.000 |
| order | 0.958 / 0.989 [0.980, 0.994] (n=1041) | 0.958 / 0.989 [0.980, 0.994] (n=1041) | +0.000 | +0.000 |
| family | 0.918 / 0.986 [0.977, 0.992] (n=1039) | 0.875 / 0.991 [0.983, 0.996] (n=1039) | +0.005 | -0.043 |
| genus | 0.874 / 0.989 [0.980, 0.994] (n=1038) | 0.872 / 0.989 [0.980, 0.994] (n=1038) | -0.000 | -0.002 |
| species | 0.727 / 1.000 [0.893, 1.000] (n=44) | 0.727 / 1.000 [0.893, 1.000] (n=44) | +0.000 | +0.000 |

**Basidiomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 1.000 / 1.000 [0.977, 1.000] (n=161) | 1.000 / 1.000 [0.977, 1.000] (n=161) | +0.000 | +0.000 |
| class | 1.000 / 0.994 [0.966, 0.999] (n=161) | 1.000 / 0.994 [0.966, 0.999] (n=161) | +0.000 | +0.000 |
| order | 1.000 / 0.975 [0.938, 0.990] (n=161) | 1.000 / 0.975 [0.938, 0.990] (n=161) | +0.000 | +0.000 |
| family | 0.969 / 0.987 [0.954, 0.996] (n=159) | 0.830 / 0.992 [0.958, 0.999] (n=159) | +0.005 | -0.138 |
| genus | 0.781 / 0.984 [0.944, 0.996] (n=160) | 0.775 / 0.984 [0.943, 0.996] (n=160) | -0.000 | -0.006 |
| species | 0.875 / 1.000 [0.646, 1.000] (n=8) | 0.875 / 1.000 [0.646, 1.000] (n=8) | +0.000 | +0.000 |

**other phyla**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.473 / 0.846 [0.665, 0.939] (n=55) | 0.473 / 0.846 [0.665, 0.939] (n=55) | +0.000 | +0.000 |
| class | 0.194 / 0.857 [0.487, 0.974] (n=36) | 0.194 / 0.857 [0.487, 0.974] (n=36) | +0.000 | +0.000 |
| order | 0.194 / 0.857 [0.487, 0.974] (n=36) | 0.194 / 0.857 [0.487, 0.974] (n=36) | +0.000 | +0.000 |
| family | 0.481 / 0.880 [0.700, 0.958] (n=52) | 0.404 / 1.000 [0.845, 1.000] (n=52) | +0.120 | -0.077 |
| genus | 0.262 / 0.909 [0.623, 0.984] (n=42) | 0.262 / 0.909 [0.623, 0.984] (n=42) | +0.000 | +0.000 |

## Truth set, leave-one-out

Cells: called share / precision [95% Wilson CI] (n with truth). The CI is on the calls made.


**all**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.975 / 0.978 [0.963, 0.987] (n=653) | 0.975 / 0.978 [0.963, 0.987] (n=653) | +0.000 | +0.000 |
| class | 0.984 / 0.971 [0.955, 0.982] (n=637) | 0.984 / 0.971 [0.955, 0.982] (n=637) | +0.000 | +0.000 |
| order | 0.973 / 0.969 [0.953, 0.980] (n=636) | 0.973 / 0.969 [0.953, 0.980] (n=636) | +0.000 | +0.000 |
| family | 0.924 / 0.986 [0.973, 0.993] (n=635) | 0.794 / 1.000 [0.992, 1.000] (n=635) | +0.014 | -0.131 |
| genus | 0.770 / 0.994 [0.983, 0.998] (n=653) | 0.767 / 0.994 [0.983, 0.998] (n=653) | -0.000 | -0.003 |
| species | 0.387 / 0.949 [0.914, 0.970] (n=653) | 0.387 / 0.949 [0.914, 0.970] (n=653) | +0.000 | +0.000 |

**Ascomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.992 / 0.988 [0.974, 0.995] (n=507) | 0.992 / 0.988 [0.974, 0.995] (n=507) | +0.000 | +0.000 |
| class | 0.990 / 0.982 [0.966, 0.991] (n=507) | 0.990 / 0.982 [0.966, 0.991] (n=507) | +0.000 | +0.000 |
| order | 0.980 / 0.982 [0.966, 0.990] (n=506) | 0.980 / 0.982 [0.966, 0.990] (n=506) | +0.000 | +0.000 |
| family | 0.944 / 0.992 [0.979, 0.997] (n=502) | 0.823 / 1.000 [0.991, 1.000] (n=502) | +0.008 | -0.122 |
| genus | 0.803 / 0.993 [0.979, 0.997] (n=507) | 0.799 / 0.993 [0.978, 0.997] (n=507) | -0.000 | -0.004 |
| species | 0.398 / 0.946 [0.905, 0.969] (n=507) | 0.398 / 0.946 [0.905, 0.969] (n=507) | +0.000 | +0.000 |

**Basidiomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.952 / 1.000 [0.963, 1.000] (n=105) | 0.952 / 1.000 [0.963, 1.000] (n=105) | +0.000 | +0.000 |
| class | 0.952 / 0.980 [0.930, 0.994] (n=105) | 0.952 / 0.980 [0.930, 0.994] (n=105) | +0.000 | +0.000 |
| order | 0.943 / 0.960 [0.901, 0.984] (n=105) | 0.943 / 0.960 [0.901, 0.984] (n=105) | +0.000 | +0.000 |
| family | 0.871 / 0.989 [0.938, 0.998] (n=101) | 0.693 / 1.000 [0.948, 1.000] (n=101) | +0.011 | -0.178 |
| genus | 0.676 / 1.000 [0.949, 1.000] (n=105) | 0.676 / 1.000 [0.949, 1.000] (n=105) | +0.000 | +0.000 |
| species | 0.314 / 0.970 [0.847, 0.995] (n=105) | 0.314 / 0.970 [0.847, 0.995] (n=105) | +0.000 | +0.000 |

**other phyla**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.829 / 0.765 [0.600, 0.876] (n=41) | 0.829 / 0.765 [0.600, 0.876] (n=41) | +0.000 | +0.000 |
| class | 1.000 / 0.720 [0.524, 0.857] (n=25) | 1.000 / 0.720 [0.524, 0.857] (n=25) | +0.000 | +0.000 |
| order | 0.960 / 0.750 [0.551, 0.880] (n=25) | 0.960 / 0.750 [0.551, 0.880] (n=25) | +0.000 | +0.000 |
| family | 0.781 / 0.880 [0.700, 0.958] (n=32) | 0.656 / 1.000 [0.845, 1.000] (n=32) | +0.120 | -0.125 |
| genus | 0.610 / 1.000 [0.867, 1.000] (n=41) | 0.610 / 1.000 [0.867, 1.000] (n=41) | +0.000 | +0.000 |
| species | 0.439 / 0.944 [0.742, 0.990] (n=41) | 0.439 / 0.944 [0.742, 0.990] (n=41) | +0.000 | +0.000 |

## Truth set, leave-species-out

Cells: called share / precision [95% Wilson CI] (n with truth). The CI is on the calls made.


**all**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.957 / 0.962 [0.943, 0.974] (n=653) | 0.957 / 0.962 [0.943, 0.974] (n=653) | +0.000 | +0.000 |
| class | 0.973 / 0.948 [0.928, 0.963] (n=637) | 0.973 / 0.948 [0.928, 0.963] (n=637) | +0.000 | +0.000 |
| order | 0.951 / 0.944 [0.922, 0.960] (n=636) | 0.951 / 0.944 [0.922, 0.960] (n=636) | +0.000 | +0.000 |
| family | 0.868 / 0.960 [0.940, 0.973] (n=635) | 0.600 / 0.995 [0.981, 0.999] (n=635) | +0.035 | -0.268 |
| genus | 0.534 / 0.980 [0.959, 0.990] (n=653) | 0.527 / 0.983 [0.962, 0.992] (n=653) | +0.003 | -0.008 |
| species | 0.064 / 0.000 [0.000, 0.084] (n=653) | 0.064 / 0.000 [0.000, 0.084] (n=653) | +0.000 | +0.000 |

**Ascomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.986 / 0.976 [0.959, 0.986] (n=507) | 0.986 / 0.976 [0.959, 0.986] (n=507) | +0.000 | +0.000 |
| class | 0.984 / 0.962 [0.941, 0.975] (n=507) | 0.984 / 0.962 [0.941, 0.975] (n=507) | +0.000 | +0.000 |
| order | 0.966 / 0.959 [0.938, 0.973] (n=506) | 0.966 / 0.959 [0.938, 0.973] (n=506) | +0.000 | +0.000 |
| family | 0.902 / 0.967 [0.946, 0.980] (n=502) | 0.657 / 0.994 [0.978, 0.998] (n=502) | +0.027 | -0.245 |
| genus | 0.594 / 0.980 [0.957, 0.991] (n=507) | 0.586 / 0.980 [0.957, 0.991] (n=507) | -0.000 | -0.008 |
| species | 0.075 / 0.000 [0.000, 0.092] (n=507) | 0.075 / 0.000 [0.000, 0.092] (n=507) | +0.000 | +0.000 |

**Basidiomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.914 / 0.990 [0.943, 0.998] (n=105) | 0.914 / 0.990 [0.943, 0.998] (n=105) | +0.000 | +0.000 |
| class | 0.914 / 0.958 [0.898, 0.984] (n=105) | 0.914 / 0.958 [0.898, 0.984] (n=105) | +0.000 | +0.000 |
| order | 0.886 / 0.935 [0.866, 0.970] (n=105) | 0.886 / 0.935 [0.866, 0.970] (n=105) | +0.000 | +0.000 |
| family | 0.782 / 0.975 [0.912, 0.993] (n=101) | 0.426 / 1.000 [0.918, 1.000] (n=101) | +0.025 | -0.356 |
| genus | 0.371 / 1.000 [0.910, 1.000] (n=105) | 0.371 / 1.000 [0.910, 1.000] (n=105) | +0.000 | +0.000 |
| species | 0.019 / 0.000 [0.000, 0.658] (n=105) | 0.019 / 0.000 [0.000, 0.658] (n=105) | +0.000 | +0.000 |

**other phyla**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.707 / 0.621 [0.440, 0.773] (n=41) | 0.707 / 0.621 [0.440, 0.773] (n=41) | +0.000 | +0.000 |
| class | 1.000 / 0.640 [0.445, 0.798] (n=25) | 1.000 / 0.640 [0.445, 0.798] (n=25) | +0.000 | +0.000 |
| order | 0.920 / 0.652 [0.449, 0.812] (n=25) | 0.920 / 0.652 [0.449, 0.812] (n=25) | +0.000 | +0.000 |
| family | 0.594 / 0.737 [0.512, 0.882] (n=32) | 0.250 / 1.000 [0.676, 1.000] (n=32) | +0.263 | -0.344 |
| genus | 0.220 / 0.889 [0.565, 0.980] (n=41) | 0.195 / 1.000 [0.676, 1.000] (n=41) | +0.111 | -0.024 |
| species | 0.049 / 0.000 [0.000, 0.658] (n=41) | 0.049 / 0.000 [0.000, 0.658] (n=41) | +0.000 | +0.000 |

## Truth set, leave-genus-out

Cells: called share / precision [95% Wilson CI] (n with truth). The CI is on the calls made.


**all**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.928 / 0.858 [0.828, 0.884] (n=653) | 0.928 / 0.858 [0.828, 0.884] (n=653) | +0.000 | +0.000 |
| class | 0.940 / 0.766 [0.731, 0.798] (n=637) | 0.940 / 0.766 [0.731, 0.798] (n=637) | +0.000 | +0.000 |
| order | 0.827 / 0.705 [0.665, 0.743] (n=636) | 0.827 / 0.705 [0.665, 0.743] (n=636) | +0.000 | +0.000 |
| family | 0.572 / 0.650 [0.600, 0.697] (n=635) | 0.117 / 0.703 [0.591, 0.795] (n=635) | +0.053 | -0.455 |
| genus | 0.072 / 0.000 [0.000, 0.076] (n=653) | 0.072 / 0.000 [0.000, 0.076] (n=653) | +0.000 | +0.000 |
| species | 0.005 / 0.000 [0.000, 0.562] (n=653) | 0.005 / 0.000 [0.000, 0.562] (n=653) | +0.000 | +0.000 |

**Ascomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.959 / 0.870 [0.838, 0.897] (n=507) | 0.959 / 0.870 [0.838, 0.897] (n=507) | +0.000 | +0.000 |
| class | 0.953 / 0.789 [0.750, 0.823] (n=507) | 0.953 / 0.789 [0.750, 0.823] (n=507) | +0.000 | +0.000 |
| order | 0.854 / 0.738 [0.695, 0.778] (n=506) | 0.854 / 0.738 [0.695, 0.778] (n=506) | +0.000 | +0.000 |
| family | 0.618 / 0.674 [0.620, 0.724] (n=502) | 0.135 / 0.721 [0.604, 0.813] (n=502) | +0.046 | -0.482 |
| genus | 0.083 / 0.000 [0.000, 0.084] (n=507) | 0.083 / 0.000 [0.000, 0.084] (n=507) | +0.000 | +0.000 |
| species | 0.006 / 0.000 [0.000, 0.562] (n=507) | 0.006 / 0.000 [0.000, 0.562] (n=507) | +0.000 | +0.000 |

**Basidiomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.876 / 0.913 [0.838, 0.955] (n=105) | 0.876 / 0.913 [0.838, 0.955] (n=105) | +0.000 | +0.000 |
| class | 0.876 / 0.728 [0.630, 0.809] (n=105) | 0.876 / 0.728 [0.630, 0.809] (n=105) | +0.000 | +0.000 |
| order | 0.714 / 0.573 [0.461, 0.679] (n=105) | 0.714 / 0.573 [0.461, 0.679] (n=105) | +0.000 | +0.000 |
| family | 0.396 / 0.525 [0.375, 0.671] (n=101) | 0.030 / 0.333 [0.061, 0.792] (n=101) | -0.192 | -0.366 |
| genus | 0.019 / 0.000 [0.000, 0.658] (n=105) | 0.019 / 0.000 [0.000, 0.658] (n=105) | +0.000 | +0.000 |
| species | 0 called (n=105) | 0 called (n=105) | +nan | +0.000 |

**other phyla**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 0.683 / 0.464 [0.295, 0.642] (n=41) | 0.683 / 0.464 [0.295, 0.642] (n=41) | +0.000 | +0.000 |
| class | 0.960 / 0.458 [0.279, 0.649] (n=25) | 0.960 / 0.458 [0.279, 0.649] (n=25) | +0.000 | +0.000 |
| order | 0.760 / 0.474 [0.273, 0.683] (n=25) | 0.760 / 0.474 [0.273, 0.683] (n=25) | +0.000 | +0.000 |
| family | 0.406 / 0.462 [0.232, 0.709] (n=32) | 0.094 / 0.667 [0.208, 0.939] (n=32) | +0.205 | -0.312 |
| genus | 0.073 / 0.000 [0.000, 0.562] (n=41) | 0.073 / 0.000 [0.000, 0.562] (n=41) | +0.000 | +0.000 |
| species | 0 called (n=41) | 0 called (n=41) | +nan | +0.000 |

## Truth set, outside the index (21)

Cells: called share / precision [95% Wilson CI] (n with truth). The CI is on the calls made.


**all**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 1.000 / 1.000 [0.845, 1.000] (n=21) | 1.000 / 1.000 [0.845, 1.000] (n=21) | +0.000 | +0.000 |
| class | 1.000 / 1.000 [0.839, 1.000] (n=20) | 1.000 / 1.000 [0.839, 1.000] (n=20) | +0.000 | +0.000 |
| order | 1.000 / 1.000 [0.839, 1.000] (n=20) | 1.000 / 1.000 [0.839, 1.000] (n=20) | +0.000 | +0.000 |
| family | 1.000 / 1.000 [0.845, 1.000] (n=21) | 1.000 / 1.000 [0.845, 1.000] (n=21) | +0.000 | +0.000 |
| genus | 1.000 / 1.000 [0.845, 1.000] (n=21) | 1.000 / 1.000 [0.845, 1.000] (n=21) | +0.000 | +0.000 |
| species | 0.857 / 1.000 [0.824, 1.000] (n=21) | 0.857 / 1.000 [0.824, 1.000] (n=21) | +0.000 | +0.000 |

**Ascomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 1.000 / 1.000 [0.816, 1.000] (n=17) | 1.000 / 1.000 [0.816, 1.000] (n=17) | +0.000 | +0.000 |
| class | 1.000 / 1.000 [0.816, 1.000] (n=17) | 1.000 / 1.000 [0.816, 1.000] (n=17) | +0.000 | +0.000 |
| order | 1.000 / 1.000 [0.816, 1.000] (n=17) | 1.000 / 1.000 [0.816, 1.000] (n=17) | +0.000 | +0.000 |
| family | 1.000 / 1.000 [0.816, 1.000] (n=17) | 1.000 / 1.000 [0.816, 1.000] (n=17) | +0.000 | +0.000 |
| genus | 1.000 / 1.000 [0.816, 1.000] (n=17) | 1.000 / 1.000 [0.816, 1.000] (n=17) | +0.000 | +0.000 |
| species | 0.824 / 1.000 [0.785, 1.000] (n=17) | 0.824 / 1.000 [0.785, 1.000] (n=17) | +0.000 | +0.000 |

**Basidiomycota**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 1.000 / 1.000 [0.438, 1.000] (n=3) | 1.000 / 1.000 [0.438, 1.000] (n=3) | +0.000 | +0.000 |
| class | 1.000 / 1.000 [0.438, 1.000] (n=3) | 1.000 / 1.000 [0.438, 1.000] (n=3) | +0.000 | +0.000 |
| order | 1.000 / 1.000 [0.438, 1.000] (n=3) | 1.000 / 1.000 [0.438, 1.000] (n=3) | +0.000 | +0.000 |
| family | 1.000 / 1.000 [0.438, 1.000] (n=3) | 1.000 / 1.000 [0.438, 1.000] (n=3) | +0.000 | +0.000 |
| genus | 1.000 / 1.000 [0.438, 1.000] (n=3) | 1.000 / 1.000 [0.438, 1.000] (n=3) | +0.000 | +0.000 |
| species | 1.000 / 1.000 [0.438, 1.000] (n=3) | 1.000 / 1.000 [0.438, 1.000] (n=3) | +0.000 | +0.000 |

**other phyla**

| rank | before | after | precision change | called change |
|---|---|---|---|---|
| phylum | 1.000 / 1.000 [0.207, 1.000] (n=1) | 1.000 / 1.000 [0.207, 1.000] (n=1) | +0.000 | +0.000 |
| family | 1.000 / 1.000 [0.207, 1.000] (n=1) | 1.000 / 1.000 [0.207, 1.000] (n=1) | +0.000 | +0.000 |
| genus | 1.000 / 1.000 [0.207, 1.000] (n=1) | 1.000 / 1.000 [0.207, 1.000] (n=1) | +0.000 | +0.000 |
| species | 1.000 / 1.000 [0.207, 1.000] (n=1) | 1.000 / 1.000 [0.207, 1.000] (n=1) | +0.000 | +0.000 |

## IndigoAg (199 assemblies, scored against provider genus names, noisy labels; index 2026.10.3)

| rank | called share before | after | precision before | after |
|---|---|---|---|---|
| family | 0.901 (154 of 171) | 0.696 (119) | 0.792 (122 correct) | 0.824 (98) |
| genus | 0.661 (115 of 174) | 0.626 (109) | 0.713 (82) | **0.706 (77)** |

## What got worse, honestly

* **Family coverage drops a lot for a small precision gain.** Post-freeze: called 0.906 to 0.850 (Basidiomycota 0.969 to 0.830), precision 0.984 to 0.992. Truth-set leave-one-out: called 0.924 to 0.794,
  precision 0.986 to 1.000. Leave-species-out: called 0.868 to 0.600. IndigoAg: called 0.90 to 0.70. On genomes with a relative in the reference the family call was already about 98-99% right, so most of the
  withheld family calls were correct ones (the 12 points of coverage the bench predicted are 6-13 points on assemblies and 27 on leave-species-out). A withheld family leaves the order call.
* **Genus: no measurable gain on real sets.** Post-freeze genus precision 0.988 to 0.988 and called 0.841 to 0.839 (3 fewer calls); leave-one-out 0.994 to 0.994 (called 0.770 to 0.767);
  leave-species-out 0.980 to 0.983 (called 0.534 to 0.527). On IndigoAg genus precision is 0.713 to 0.706, 5 fewer correct calls and 6 fewer calls (not significant, but not an improvement).
* Basidiomycota genus is not worse after the change on the held-out bench (0.955 to 0.960) and unchanged post-freeze (0.984), but the Wilson interval there is wide (0.944-0.996).
* Leave-genus-out family for Basidiomycota: 0.525 (40 calls) to 0.333 (3 calls), and for all phyla 0.650 to 0.703 with the called share falling from 0.572 to 0.117; those are the inputs for which the genus is absent,
  where a family name is the most that can be said, so the call share is the cost.
* Held-out bench coverage in the derivation is in-sample for the family cutoff (fit on the full pool); the 10-split numbers above are out-of-sample.

## Recommendation

The family margin does what it was designed to do (precision up, calls down) but costs more calls on real assemblies than the bench suggested. The genus change is close to a no-op because #57 already
raised the genus cutoff above what a 0.97 fit asks for. If the family coverage loss is not wanted, revert only the family row (`min_ani[3]` back to 0.7715): genus is unaffected by it.
