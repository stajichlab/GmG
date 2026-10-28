# Clade-aware calibration strategies (target 0.95)

70878 bench queries, 10 random genus splits; fit on half the genera, score on the other half, pooled over splits. Pool `known` (leave-one-out for genus/species, leave-species-out above). `phylum` fits a threshold per predicted phylum when it has ≥ 400 voted queries.


## True phylum group: all

| rank | global prec (cov) | margin prec (cov) | weighted prec (cov) | phylum prec (cov) |
|---|---|---|---|---|
| phylum | 0.961 (1.00) | 0.965 (0.99) | 0.998 (0.75) | 0.962 (0.99) |
| class | 0.949 (0.99) | 0.959 (0.97) | 0.986 (0.80) | 0.950 (0.99) |
| order | 0.942 (0.98) | 0.957 (0.93) | 0.978 (0.78) | 0.943 (0.98) |
| family | 0.929 (0.96) | 0.946 (0.82) | 0.954 (0.33) | 0.926 (0.97) |
| genus | 0.944 (0.98) | 0.965 (0.94) | 0.987 (0.74) | 0.947 (0.97) |
| species | 0.950 (0.55) | 0.948 (0.11) | 0.957 (0.32) | 0.947 (0.54) |

## True phylum group: Ascomycota

| rank | global prec (cov) | margin prec (cov) | weighted prec (cov) | phylum prec (cov) |
|---|---|---|---|---|
| phylum | 0.970 (1.00) | 0.973 (0.99) | 0.999 (0.77) | 0.971 (1.00) |
| class | 0.952 (0.99) | 0.959 (0.98) | 0.984 (0.82) | 0.953 (0.99) |
| order | 0.945 (0.99) | 0.956 (0.94) | 0.976 (0.79) | 0.944 (0.99) |
| family | 0.930 (0.97) | 0.944 (0.84) | 0.951 (0.35) | 0.924 (0.98) |
| genus | 0.956 (0.99) | 0.972 (0.96) | 0.990 (0.78) | 0.951 (0.99) |
| species | 0.957 (0.61) | 0.952 (0.11) | 0.962 (0.36) | 0.953 (0.63) |

## True phylum group: Basidiomycota

| rank | global prec (cov) | margin prec (cov) | weighted prec (cov) | phylum prec (cov) |
|---|---|---|---|---|
| phylum | 0.961 (0.99) | 0.968 (0.98) | 0.999 (0.65) | 0.960 (0.95) |
| class | 0.949 (0.98) | 0.968 (0.95) | 0.997 (0.73) | 0.951 (0.95) |
| order | 0.938 (0.96) | 0.970 (0.89) | 0.986 (0.69) | 0.944 (0.94) |
| family | 0.947 (0.93) | 0.974 (0.75) | 0.979 (0.18) | 0.953 (0.89) |
| genus | 0.910 (0.96) | 0.937 (0.89) | 0.969 (0.54) | 0.943 (0.84) |
| species | 0.926 (0.33) | 0.947 (0.10) | 0.936 (0.22) | 0.925 (0.14) |

## True phylum group: other

| rank | global prec (cov) | margin prec (cov) | weighted prec (cov) | phylum prec (cov) |
|---|---|---|---|---|
| phylum | 0.809 (0.99) | 0.835 (0.96) | 0.967 (0.71) | 0.819 (0.98) |
| class | 0.884 (0.98) | 0.914 (0.94) | 0.986 (0.80) | 0.893 (0.97) |
| order | 0.898 (0.97) | 0.939 (0.91) | 0.987 (0.79) | 0.906 (0.96) |
| family | 0.867 (0.92) | 0.905 (0.79) | 0.970 (0.34) | 0.874 (0.91) |
| genus | 0.860 (0.96) | 0.915 (0.88) | 0.968 (0.59) | 0.878 (0.94) |
| species | 0.817 (0.34) | 0.859 (0.08) | 0.876 (0.12) | 0.817 (0.34) |

## Reading

* The shipped approach (`global`) delivers 0.956 genus precision on held-out Ascomycota genera but 0.910 on Basidiomycota and 0.860 on the remaining phyla, against a 0.95 target.
* `weighted` (equal weight per phylum) meets the target everywhere but fits to the worst clade, so coverage drops to 0.74 overall (family 0.33): too costly as a default.
* `margin` (fit at 0.97) gains 2-5 points in the weaker clades for 4-6 points of coverage, but collapses species coverage (0.55 → 0.11), because the species cutoff sits on a cliff.
* `phylum` (per predicted phylum, ≥ 400 queries) helps where there is enough data: Basidiomycota genus 0.910 → 0.943 (coverage 0.96 → 0.84), and does nothing for the small phyla, which fall back to the global fit (0.86 genus, 0.81 phylum).
* No strategy fixes the small phyla: they have too few genomes to calibrate, and their errors are mostly neighbours in other phyla. This is a reference-coverage limit, not a threshold one.
* These runs fit the `known` pool at the old genus cutoff (~0.755); the shipped balanced default (genus ≥ 0.90) already raises genus precision, so the Basidiomycota gap is smaller there (see `post_freeze_balanced.md`).
