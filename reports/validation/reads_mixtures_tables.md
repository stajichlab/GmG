# Mixed read sets (#45): all scored tables

Output of `scripts/mix_score.py` for the runs of `scripts/hpcc/34_mixtures.sbatch` (default read settings and `--min-abund 1`; 517 inputs each). The summary and the reading of these tables are in [reads_mixtures.md](reads_mixtures.md).

## Run `default`: 513 two-fungus inputs


### Depth `high`

| distance | rank | A share | n | majority / shared | minority | other | no call | minority in top-10 hits | ANI vs A alone: observed / p^(1/21) | genus of A kept | weak_phylum |
|---|---|---|---|---|---|---|---|---|---|---|---|
| same_species | genus | 100% | 6 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_species | genus | 99.9% | 6 | 100% | 0% | 0% | 0% | - | 0.9999 / 1.0000 | 100% | 0 |
| same_species | genus | 99.5% | 6 | 100% | 0% | 0% | 0% | - | 0.9997 / 0.9998 | 100% | 0 |
| same_species | genus | 99% | 6 | 100% | 0% | 0% | 0% | - | 0.9994 / 0.9995 | 100% | 0 |
| same_species | genus | 95% | 6 | 100% | 0% | 0% | 0% | - | 0.9964 / 0.9976 | 100% | 0 |
| same_species | genus | 90% | 6 | 100% | 0% | 0% | 0% | - | 0.9936 / 0.9950 | 100% | 0 |
| same_species | genus | 80% | 6 | 100% | 0% | 0% | 0% | - | 0.9913 / 0.9894 | 100% | 0 |
| same_species | genus | 50% | 6 | 100% | 0% | 0% | 0% | - | 0.9887 / 0.9675 | 100% | 0 |
| same_species | genus | 0% | 6 | 100% | 0% | 0% | 0% | - | - | 100% | 0 |
| same_genus | genus | 100% | 7 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_genus | genus | 99.9% | 7 | 100% | 0% | 0% | 0% | - | 0.9999 / 1.0000 | 100% | 0 |
| same_genus | genus | 99.5% | 7 | 100% | 0% | 0% | 0% | - | 0.9996 / 0.9998 | 100% | 0 |
| same_genus | genus | 99% | 7 | 100% | 0% | 0% | 0% | - | 0.9992 / 0.9995 | 100% | 0 |
| same_genus | genus | 95% | 7 | 100% | 0% | 0% | 0% | - | 0.9999 / 0.9976 | 100% | 0 |
| same_genus | genus | 90% | 7 | 100% | 0% | 0% | 0% | - | 1.0004 / 0.9950 | 100% | 0 |
| same_genus | genus | 80% | 7 | 100% | 0% | 0% | 0% | - | 0.9998 / 0.9894 | 100% | 0 |
| same_genus | genus | 50% | 7 | 100% | 0% | 0% | 0% | - | 0.9960 / 0.9675 | 100% | 0 |
| same_genus | genus | 0% | 7 | 100% | 0% | 0% | 0% | - | - | 100% | 0 |
| same_order | family | 100% | 7 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_order | family | 99.9% | 7 | 100% | 0% | 0% | 0% | 0% | 0.9999 / 1.0000 | 100% | 0 |
| same_order | family | 99.5% | 7 | 100% | 0% | 0% | 0% | 40% | 0.9992 / 0.9998 | 100% | 0 |
| same_order | family | 99% | 7 | 100% | 0% | 0% | 0% | 40% | 0.9984 / 0.9995 | 100% | 0 |
| same_order | family | 95% | 7 | 100% | 0% | 0% | 0% | 80% | 0.9929 / 0.9976 | 100% | 0 |
| same_order | family | 90% | 7 | 71% | 0% | 0% | 29% | 80% | 0.9916 / 0.9950 | 57% | 0 |
| same_order | family | 80% | 7 | 43% | 29% | 0% | 29% | 100% | 0.9955 / 0.9894 | 29% | 0 |
| same_order | family | 50% | 7 | 100% | 0% | 0% | 0% | 100% | 1.0045 / 0.9675 | 29% | 0 |
| same_order | family | 0% | 7 | 100% | 0% | 0% | 0% | - | - | 0% | 0 |
| diff_phylum | phylum | 100% | 5 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| diff_phylum | phylum | 99.9% | 5 | 100% | 0% | 0% | 0% | 20% | 0.9994 / 1.0000 | 100% | 0 |
| diff_phylum | phylum | 99.5% | 5 | 100% | 0% | 0% | 0% | 100% | 0.9975 / 0.9998 | 100% | 0 |
| diff_phylum | phylum | 99% | 5 | 100% | 0% | 0% | 0% | 100% | 0.9958 / 0.9995 | 100% | 0 |
| diff_phylum | phylum | 95% | 5 | 80% | 20% | 0% | 0% | 100% | 0.9933 / 0.9976 | 80% | 0 |
| diff_phylum | phylum | 90% | 5 | 80% | 20% | 0% | 0% | 100% | 0.9912 / 0.9950 | 80% | 0 |
| diff_phylum | phylum | 80% | 5 | 80% | 20% | 0% | 0% | 100% | 0.9862 / 0.9894 | 80% | 0 |
| diff_phylum | phylum | 50% | 5 | 100% | 0% | 0% | 0% | 100% | 0.9816 / 0.9675 | 40% | 0 |
| diff_phylum | phylum | 0% | 5 | 100% | 0% | 0% | 0% | - | - | 0% | 0 |

### Depth `low`

| distance | rank | A share | n | majority / shared | minority | other | no call | minority in top-10 hits | ANI vs A alone: observed / p^(1/21) | genus of A kept | weak_phylum |
|---|---|---|---|---|---|---|---|---|---|---|---|
| same_species | genus | 100% | 8 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_species | genus | 99.9% | 8 | 100% | 0% | 0% | 0% | - | 1.0001 / 1.0000 | 100% | 0 |
| same_species | genus | 99.5% | 8 | 100% | 0% | 0% | 0% | - | 0.9999 / 0.9998 | 100% | 0 |
| same_species | genus | 99% | 8 | 100% | 0% | 0% | 0% | - | 1.0000 / 0.9995 | 100% | 0 |
| same_species | genus | 95% | 8 | 100% | 0% | 0% | 0% | - | 0.9998 / 0.9976 | 100% | 0 |
| same_species | genus | 90% | 8 | 100% | 0% | 0% | 0% | - | 0.9993 / 0.9950 | 100% | 0 |
| same_species | genus | 80% | 8 | 100% | 0% | 0% | 0% | - | 0.9977 / 0.9894 | 100% | 0 |
| same_species | genus | 50% | 8 | 100% | 0% | 0% | 0% | - | 0.9934 / 0.9675 | 100% | 0 |
| same_species | genus | 0% | 8 | 100% | 0% | 0% | 0% | - | - | 100% | 0 |
| same_genus | genus | 100% | 8 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_genus | genus | 99.9% | 8 | 100% | 0% | 0% | 0% | - | 0.9999 / 1.0000 | 100% | 0 |
| same_genus | genus | 99.5% | 8 | 100% | 0% | 0% | 0% | - | 0.9997 / 0.9998 | 100% | 0 |
| same_genus | genus | 99% | 8 | 100% | 0% | 0% | 0% | - | 0.9994 / 0.9995 | 100% | 0 |
| same_genus | genus | 95% | 8 | 100% | 0% | 0% | 0% | - | 0.9981 / 0.9976 | 100% | 0 |
| same_genus | genus | 90% | 8 | 100% | 0% | 0% | 0% | - | 0.9987 / 0.9950 | 100% | 0 |
| same_genus | genus | 80% | 8 | 100% | 0% | 0% | 0% | - | 0.9998 / 0.9894 | 100% | 0 |
| same_genus | genus | 50% | 8 | 100% | 0% | 0% | 0% | - | 0.9926 / 0.9675 | 100% | 0 |
| same_genus | genus | 0% | 8 | 100% | 0% | 0% | 0% | - | - | 100% | 0 |
| same_order | family | 100% | 8 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_order | family | 99.9% | 8 | 100% | 0% | 0% | 0% | 17% | 0.9999 / 1.0000 | 100% | 0 |
| same_order | family | 99.5% | 8 | 100% | 0% | 0% | 0% | 33% | 0.9996 / 0.9998 | 100% | 0 |
| same_order | family | 99% | 8 | 100% | 0% | 0% | 0% | 50% | 0.9990 / 0.9995 | 100% | 0 |
| same_order | family | 95% | 8 | 100% | 0% | 0% | 0% | 67% | 0.9993 / 0.9976 | 88% | 0 |
| same_order | family | 90% | 8 | 100% | 0% | 0% | 0% | 83% | 0.9993 / 0.9950 | 88% | 0 |
| same_order | family | 80% | 8 | 75% | 12% | 0% | 12% | 83% | 0.9996 / 0.9894 | 62% | 0 |
| same_order | family | 50% | 8 | 100% | 0% | 0% | 0% | 100% | 1.0096 / 0.9675 | 38% | 0 |
| same_order | family | 0% | 8 | 100% | 0% | 0% | 0% | - | - | 0% | 0 |
| diff_phylum | phylum | 100% | 8 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| diff_phylum | phylum | 99.9% | 8 | 100% | 0% | 0% | 0% | 12% | 0.9996 / 1.0000 | 100% | 0 |
| diff_phylum | phylum | 99.5% | 8 | 100% | 0% | 0% | 0% | 62% | 0.9985 / 0.9998 | 100% | 0 |
| diff_phylum | phylum | 99% | 8 | 100% | 0% | 0% | 0% | 75% | 0.9971 / 0.9995 | 100% | 0 |
| diff_phylum | phylum | 95% | 8 | 88% | 12% | 0% | 0% | 88% | 0.9955 / 0.9976 | 88% | 0 |
| diff_phylum | phylum | 90% | 8 | 88% | 12% | 0% | 0% | 100% | 0.9937 / 0.9950 | 88% | 0 |
| diff_phylum | phylum | 80% | 8 | 75% | 25% | 0% | 0% | 100% | 0.9889 / 0.9894 | 62% | 0 |
| diff_phylum | phylum | 50% | 8 | 100% | 0% | 0% | 0% | 100% | 0.9937 / 0.9675 | 38% | 0 |
| diff_phylum | phylum | 0% | 8 | 100% | 0% | 0% | 0% | - | - | 0% | 0 |

**Coccidioides pairs (curated species): the species called**

| distance | depth | A share | n | species of A | species of B | other or none |
|---|---|---|---|---|---|---|
| same_species | high | 100% | 2 | 0% | 0% | 100% |
| same_species | high | 99.9% | 2 | 0% | 0% | 100% |
| same_species | high | 99.5% | 2 | 0% | 0% | 100% |
| same_species | high | 99% | 2 | 0% | 0% | 100% |
| same_species | high | 95% | 2 | 0% | 0% | 100% |
| same_species | high | 90% | 2 | 0% | 0% | 100% |
| same_species | high | 80% | 2 | 0% | 0% | 100% |
| same_species | high | 50% | 2 | 0% | 0% | 100% |
| same_species | high | 0% | 2 | 50% | 0% | 50% |
| same_species | low | 100% | 2 | 50% | 0% | 50% |
| same_species | low | 99.9% | 2 | 50% | 0% | 50% |
| same_species | low | 99.5% | 2 | 50% | 0% | 50% |
| same_species | low | 99% | 2 | 50% | 0% | 50% |
| same_species | low | 95% | 2 | 50% | 0% | 50% |
| same_species | low | 90% | 2 | 50% | 0% | 50% |
| same_species | low | 80% | 2 | 50% | 0% | 50% |
| same_species | low | 50% | 2 | 50% | 0% | 50% |
| same_species | low | 0% | 2 | 100% | 0% | 0% |
| same_genus | high | 100% | 1 | 0% | 0% | 100% |
| same_genus | high | 99.9% | 1 | 0% | 0% | 100% |
| same_genus | high | 99.5% | 1 | 0% | 0% | 100% |
| same_genus | high | 99% | 1 | 0% | 0% | 100% |
| same_genus | high | 95% | 1 | 0% | 0% | 100% |
| same_genus | high | 90% | 1 | 0% | 0% | 100% |
| same_genus | high | 80% | 1 | 0% | 0% | 100% |
| same_genus | high | 50% | 1 | 0% | 0% | 100% |
| same_genus | high | 0% | 1 | 0% | 0% | 100% |
| same_genus | low | 100% | 1 | 0% | 0% | 100% |
| same_genus | low | 99.9% | 1 | 0% | 0% | 100% |
| same_genus | low | 99.5% | 1 | 0% | 0% | 100% |
| same_genus | low | 99% | 1 | 0% | 0% | 100% |
| same_genus | low | 95% | 1 | 0% | 0% | 100% |
| same_genus | low | 90% | 1 | 0% | 0% | 100% |
| same_genus | low | 80% | 1 | 0% | 0% | 100% |
| same_genus | low | 50% | 1 | 0% | 0% | 100% |
| same_genus | low | 0% | 1 | 0% | 100% | 0% |

**Three-genome mixtures (60:30:10)**

| triple | type | A (60%) | B (30%) | C (10%) | deepest call | genus call is A / B / C | genus of A / B / C in the top-10 hits |
|---|---|---|---|---|---|---|---|
| t00 | same_order | Lichtheimia | Backusella | Apophysomyces | Lichtheimia | yes/no/no | yes/yes/no |
| t01 | same_order | Lichtheimia | Benjaminiella | Pilaira | Lichtheimia | yes/no/no | yes/yes/yes |
| t02 | three_phyla | Coccidioides | Phlyctochytrium | Kickxella | Coccidioides | yes/no/no | yes/no/no |
| t03 | three_phyla | Coccidioides | Entophlyctis | Dispira | Coccidioides | yes/no/no | yes/no/no |

## Run `abund1`: 513 two-fungus inputs


### Depth `high`

| distance | rank | A share | n | majority / shared | minority | other | no call | minority in top-10 hits | ANI vs A alone: observed / p^(1/21) | genus of A kept | weak_phylum |
|---|---|---|---|---|---|---|---|---|---|---|---|
| same_species | genus | 100% | 6 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_species | genus | 99.9% | 6 | 100% | 0% | 0% | 0% | - | 0.9998 / 1.0000 | 100% | 0 |
| same_species | genus | 99.5% | 6 | 100% | 0% | 0% | 0% | - | 0.9994 / 0.9998 | 100% | 0 |
| same_species | genus | 99% | 6 | 100% | 0% | 0% | 0% | - | 0.9988 / 0.9995 | 100% | 0 |
| same_species | genus | 95% | 6 | 100% | 0% | 0% | 0% | - | 0.9965 / 0.9976 | 100% | 0 |
| same_species | genus | 90% | 6 | 100% | 0% | 0% | 0% | - | 0.9956 / 0.9950 | 100% | 0 |
| same_species | genus | 80% | 6 | 100% | 0% | 0% | 0% | - | 0.9952 / 0.9894 | 100% | 0 |
| same_species | genus | 50% | 6 | 100% | 0% | 0% | 0% | - | 0.9952 / 0.9675 | 100% | 0 |
| same_species | genus | 0% | 6 | 100% | 0% | 0% | 0% | - | - | 100% | 0 |
| same_genus | genus | 100% | 7 | 86% | 0% | 0% | 14% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_genus | genus | 99.9% | 7 | 86% | 0% | 0% | 14% | - | 0.9998 / 1.0000 | 100% | 0 |
| same_genus | genus | 99.5% | 7 | 86% | 0% | 0% | 14% | - | 0.9995 / 0.9998 | 100% | 0 |
| same_genus | genus | 99% | 7 | 86% | 0% | 0% | 14% | - | 0.9991 / 0.9995 | 100% | 0 |
| same_genus | genus | 95% | 7 | 100% | 0% | 0% | 0% | - | 1.0027 / 0.9976 | 100% | 0 |
| same_genus | genus | 90% | 7 | 100% | 0% | 0% | 0% | - | 1.0044 / 0.9950 | 100% | 0 |
| same_genus | genus | 80% | 7 | 100% | 0% | 0% | 0% | - | 1.0054 / 0.9894 | 100% | 0 |
| same_genus | genus | 50% | 7 | 100% | 0% | 0% | 0% | - | 1.0029 / 0.9675 | 100% | 0 |
| same_genus | genus | 0% | 7 | 100% | 0% | 0% | 0% | - | - | 100% | 0 |
| same_order | family | 100% | 7 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_order | family | 99.9% | 7 | 100% | 0% | 0% | 0% | 0% | 0.9998 / 1.0000 | 100% | 0 |
| same_order | family | 99.5% | 7 | 100% | 0% | 0% | 0% | 40% | 0.9991 / 0.9998 | 100% | 0 |
| same_order | family | 99% | 7 | 100% | 0% | 0% | 0% | 40% | 0.9983 / 0.9995 | 100% | 0 |
| same_order | family | 95% | 7 | 100% | 0% | 0% | 0% | 80% | 0.9944 / 0.9976 | 83% | 0 |
| same_order | family | 90% | 7 | 71% | 29% | 0% | 0% | 80% | 0.9976 / 0.9950 | 67% | 0 |
| same_order | family | 80% | 7 | 43% | 57% | 0% | 0% | 100% | 1.0029 / 0.9894 | 33% | 0 |
| same_order | family | 50% | 7 | 100% | 0% | 0% | 0% | 100% | 1.0066 / 0.9675 | 33% | 0 |
| same_order | family | 0% | 7 | 100% | 0% | 0% | 0% | - | - | 0% | 0 |
| diff_phylum | phylum | 100% | 5 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| diff_phylum | phylum | 99.9% | 5 | 100% | 0% | 0% | 0% | 40% | 0.9996 / 1.0000 | 100% | 0 |
| diff_phylum | phylum | 99.5% | 5 | 100% | 0% | 0% | 0% | 80% | 0.9983 / 0.9998 | 100% | 0 |
| diff_phylum | phylum | 99% | 5 | 100% | 0% | 0% | 0% | 100% | 0.9970 / 0.9995 | 100% | 0 |
| diff_phylum | phylum | 95% | 5 | 80% | 20% | 0% | 0% | 100% | 0.9953 / 0.9976 | 75% | 0 |
| diff_phylum | phylum | 90% | 5 | 80% | 20% | 0% | 0% | 100% | 0.9946 / 0.9950 | 75% | 0 |
| diff_phylum | phylum | 80% | 5 | 60% | 40% | 0% | 0% | 100% | 0.9934 / 0.9894 | 50% | 0 |
| diff_phylum | phylum | 50% | 5 | 100% | 0% | 0% | 0% | 100% | 0.9927 / 0.9675 | 25% | 0 |
| diff_phylum | phylum | 0% | 5 | 100% | 0% | 0% | 0% | - | - | 0% | 0 |

### Depth `low`

| distance | rank | A share | n | majority / shared | minority | other | no call | minority in top-10 hits | ANI vs A alone: observed / p^(1/21) | genus of A kept | weak_phylum |
|---|---|---|---|---|---|---|---|---|---|---|---|
| same_species | genus | 100% | 8 | 88% | 0% | 0% | 12% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_species | genus | 99.9% | 8 | 88% | 0% | 0% | 12% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_species | genus | 99.5% | 8 | 88% | 0% | 0% | 12% | - | 0.9997 / 0.9998 | 100% | 0 |
| same_species | genus | 99% | 8 | 88% | 0% | 0% | 12% | - | 0.9995 / 0.9995 | 100% | 0 |
| same_species | genus | 95% | 8 | 88% | 0% | 0% | 12% | - | 0.9989 / 0.9976 | 100% | 0 |
| same_species | genus | 90% | 8 | 100% | 0% | 0% | 0% | - | 0.9983 / 0.9950 | 100% | 0 |
| same_species | genus | 80% | 8 | 100% | 0% | 0% | 0% | - | 0.9973 / 0.9894 | 100% | 0 |
| same_species | genus | 50% | 8 | 100% | 0% | 0% | 0% | - | 0.9958 / 0.9675 | 100% | 0 |
| same_species | genus | 0% | 8 | 100% | 0% | 0% | 0% | - | - | 100% | 0 |
| same_genus | genus | 100% | 8 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_genus | genus | 99.9% | 8 | 88% | 0% | 0% | 12% | - | 0.9999 / 1.0000 | 88% | 0 |
| same_genus | genus | 99.5% | 8 | 88% | 0% | 0% | 12% | - | 0.9996 / 0.9998 | 88% | 0 |
| same_genus | genus | 99% | 8 | 88% | 0% | 0% | 12% | - | 0.9992 / 0.9995 | 88% | 0 |
| same_genus | genus | 95% | 8 | 100% | 0% | 0% | 0% | - | 0.9993 / 0.9976 | 100% | 0 |
| same_genus | genus | 90% | 8 | 100% | 0% | 0% | 0% | - | 1.0006 / 0.9950 | 100% | 0 |
| same_genus | genus | 80% | 8 | 100% | 0% | 0% | 0% | - | 1.0009 / 0.9894 | 100% | 0 |
| same_genus | genus | 50% | 8 | 75% | 0% | 0% | 25% | - | 0.9945 / 0.9675 | 75% | 0 |
| same_genus | genus | 0% | 8 | 88% | 0% | 0% | 12% | - | - | 88% | 0 |
| same_order | family | 100% | 8 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| same_order | family | 99.9% | 8 | 100% | 0% | 0% | 0% | 17% | 0.9998 / 1.0000 | 100% | 0 |
| same_order | family | 99.5% | 8 | 100% | 0% | 0% | 0% | 33% | 0.9993 / 0.9998 | 100% | 0 |
| same_order | family | 99% | 8 | 100% | 0% | 0% | 0% | 50% | 0.9987 / 0.9995 | 100% | 0 |
| same_order | family | 95% | 8 | 100% | 0% | 0% | 0% | 67% | 0.9984 / 0.9976 | 88% | 0 |
| same_order | family | 90% | 8 | 88% | 0% | 0% | 12% | 83% | 0.9981 / 0.9950 | 75% | 0 |
| same_order | family | 80% | 8 | 75% | 25% | 0% | 0% | 83% | 0.9997 / 0.9894 | 62% | 0 |
| same_order | family | 50% | 8 | 100% | 0% | 0% | 0% | 100% | 1.0077 / 0.9675 | 38% | 0 |
| same_order | family | 0% | 8 | 100% | 0% | 0% | 0% | - | - | 0% | 0 |
| diff_phylum | phylum | 100% | 8 | 100% | 0% | 0% | 0% | - | 1.0000 / 1.0000 | 100% | 0 |
| diff_phylum | phylum | 99.9% | 8 | 100% | 0% | 0% | 0% | 12% | 0.9996 / 1.0000 | 88% | 0 |
| diff_phylum | phylum | 99.5% | 8 | 100% | 0% | 0% | 0% | 75% | 0.9983 / 0.9998 | 88% | 0 |
| diff_phylum | phylum | 99% | 8 | 100% | 0% | 0% | 0% | 88% | 0.9968 / 0.9995 | 88% | 0 |
| diff_phylum | phylum | 95% | 8 | 88% | 12% | 0% | 0% | 100% | 0.9922 / 0.9976 | 62% | 0 |
| diff_phylum | phylum | 90% | 8 | 50% | 50% | 0% | 0% | 100% | 0.9946 / 0.9950 | 38% | 0 |
| diff_phylum | phylum | 80% | 8 | 50% | 50% | 0% | 0% | 100% | 0.9996 / 0.9894 | 25% | 0 |
| diff_phylum | phylum | 50% | 8 | 100% | 0% | 0% | 0% | 100% | 1.0054 / 0.9675 | 25% | 0 |
| diff_phylum | phylum | 0% | 8 | 100% | 0% | 0% | 0% | - | - | 0% | 0 |

**Coccidioides pairs (curated species): the species called**

| distance | depth | A share | n | species of A | species of B | other or none |
|---|---|---|---|---|---|---|
| same_species | high | 100% | 2 | 0% | 0% | 100% |
| same_species | high | 99.9% | 2 | 0% | 0% | 100% |
| same_species | high | 99.5% | 2 | 0% | 0% | 100% |
| same_species | high | 99% | 2 | 0% | 0% | 100% |
| same_species | high | 95% | 2 | 0% | 0% | 100% |
| same_species | high | 90% | 2 | 0% | 0% | 100% |
| same_species | high | 80% | 2 | 0% | 0% | 100% |
| same_species | high | 50% | 2 | 0% | 0% | 100% |
| same_species | high | 0% | 2 | 0% | 0% | 100% |
| same_species | low | 100% | 2 | 50% | 0% | 50% |
| same_species | low | 99.9% | 2 | 50% | 0% | 50% |
| same_species | low | 99.5% | 2 | 50% | 0% | 50% |
| same_species | low | 99% | 2 | 50% | 0% | 50% |
| same_species | low | 95% | 2 | 50% | 0% | 50% |
| same_species | low | 90% | 2 | 50% | 0% | 50% |
| same_species | low | 80% | 2 | 50% | 0% | 50% |
| same_species | low | 50% | 2 | 50% | 0% | 50% |
| same_species | low | 0% | 2 | 50% | 0% | 50% |
| same_genus | high | 100% | 1 | 0% | 0% | 100% |
| same_genus | high | 99.9% | 1 | 0% | 0% | 100% |
| same_genus | high | 99.5% | 1 | 0% | 0% | 100% |
| same_genus | high | 99% | 1 | 0% | 0% | 100% |
| same_genus | high | 95% | 1 | 0% | 0% | 100% |
| same_genus | high | 90% | 1 | 0% | 0% | 100% |
| same_genus | high | 80% | 1 | 0% | 0% | 100% |
| same_genus | high | 50% | 1 | 0% | 0% | 100% |
| same_genus | high | 0% | 1 | 0% | 0% | 100% |
| same_genus | low | 100% | 1 | 0% | 0% | 100% |
| same_genus | low | 99.9% | 1 | 0% | 0% | 100% |
| same_genus | low | 99.5% | 1 | 0% | 0% | 100% |
| same_genus | low | 99% | 1 | 0% | 0% | 100% |
| same_genus | low | 95% | 1 | 0% | 0% | 100% |
| same_genus | low | 90% | 1 | 0% | 0% | 100% |
| same_genus | low | 80% | 1 | 0% | 0% | 100% |
| same_genus | low | 50% | 1 | 0% | 0% | 100% |
| same_genus | low | 0% | 1 | 0% | 100% | 0% |

**Three-genome mixtures (60:30:10)**

| triple | type | A (60%) | B (30%) | C (10%) | deepest call | genus call is A / B / C | genus of A / B / C in the top-10 hits |
|---|---|---|---|---|---|---|---|
| t00 | same_order | Lichtheimia | Backusella | Apophysomyces | Backusella | no/yes/no | yes/yes/no |
| t01 | same_order | Lichtheimia | Benjaminiella | Pilaira | Lichtheimia | yes/no/no | yes/yes/yes |
| t02 | three_phyla | Coccidioides | Phlyctochytrium | Kickxella | Onygenales | no/no/no | yes/yes/no |
| t03 | three_phyla | Coccidioides | Entophlyctis | Dispira | Coccidioides | yes/no/no | yes/no/no |
