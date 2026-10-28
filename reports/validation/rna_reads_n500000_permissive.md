# Real RNA-seq reads (BFD, normalised, species-labelled): first 500000 reads of R1; DNA genome index, earlier permissive thresholds

345 inputs scored (0 without a truth row); status: {'ok': np.int64(339), 'low_confidence': np.int64(5), 'no_match': np.int64(1)}

Precision = correct among the calls made; called = share of inputs with a call at that rank; recall = correct calls / inputs with truth. Errors are excluded. Truth = the lineage in the truth table (NCBI taxonomy or curated); where it is blank the rank is not scored.

## All inputs

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 345 | 0.997 | 0.988 [0.970, 0.995] | 0.986 |
| class | 336 | 1.000 | 0.982 [0.962, 0.992] | 0.982 |
| order | 334 | 0.997 | 0.985 [0.965, 0.994] | 0.982 |
| family | 335 | 0.988 | 0.988 [0.969, 0.995] | 0.976 |
| genus | 342 | 0.985 | 0.988 [0.970, 0.995] | 0.974 |
| species | 345 | 0.594 | 0.971 [0.938, 0.987] | 0.577 |

## By `phylum`


**phylum = Ascomycota** (n = 188)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 188 | 1.000 | 0.984 [0.954, 0.995] | 0.984 |
| class | 188 | 1.000 | 0.979 [0.947, 0.992] | 0.979 |
| order | 187 | 1.000 | 0.979 [0.946, 0.992] | 0.979 |
| family | 183 | 0.995 | 0.978 [0.945, 0.991] | 0.973 |
| genus | 185 | 0.995 | 0.978 [0.945, 0.992] | 0.973 |
| species | 188 | 0.702 | 0.977 [0.935, 0.992] | 0.686 |

**phylum = Basidiomycota** (n = 91)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 91 | 1.000 | 0.989 [0.940, 0.998] | 0.989 |
| class | 91 | 1.000 | 0.978 [0.923, 0.994] | 0.978 |
| order | 91 | 0.989 | 0.989 [0.940, 0.998] | 0.978 |
| family | 89 | 0.978 | 1.000 [0.958, 1.000] | 0.978 |
| genus | 91 | 0.978 | 1.000 [0.959, 1.000] | 0.978 |
| species | 91 | 0.473 | 0.953 [0.845, 0.987] | 0.451 |

**phylum = Chytridiomycota** (n = 8)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| class | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| order | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| family | 7 | 1.000 | 1.000 [0.646, 1.000] | 1.000 |
| genus | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| species | 8 | 0.750 | 0.833 [0.436, 0.970] | 0.625 |

**phylum = Glomeromycota** (n = 8)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| class | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| order | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| family | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| genus | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| species | 8 | 0.250 | 1.000 [0.342, 1.000] | 0.250 |

**phylum = Kickxellomycota** (n = 8)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| class | 7 | 1.000 | 1.000 [0.646, 1.000] | 1.000 |
| order | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| family | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| genus | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| species | 8 | 0.625 | 1.000 [0.566, 1.000] | 0.625 |

**phylum = Microsporidia** (n = 8)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 8 | 0.875 | 1.000 [0.646, 1.000] | 0.875 |
| family | 7 | 1.000 | 1.000 [0.646, 1.000] | 1.000 |
| genus | 8 | 0.875 | 1.000 [0.646, 1.000] | 0.875 |
| species | 8 | 0.000 | - | 0.000 |

**phylum = Mortierellomycota** (n = 8)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| class | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| order | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| family | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| genus | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| species | 8 | 0.500 | 1.000 [0.510, 1.000] | 0.500 |

**phylum = Mucoromycota** (n = 8)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| class | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| order | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| family | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| genus | 8 | 1.000 | 1.000 [0.676, 1.000] | 1.000 |
| species | 8 | 0.875 | 1.000 [0.646, 1.000] | 0.875 |

**phylum = Neocallimastigomycota** (n = 6)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 6 | 1.000 | 1.000 [0.610, 1.000] | 1.000 |
| class | 6 | 1.000 | 1.000 [0.610, 1.000] | 1.000 |
| order | 6 | 1.000 | 1.000 [0.610, 1.000] | 1.000 |
| family | 6 | 1.000 | 1.000 [0.610, 1.000] | 1.000 |
| genus | 6 | 1.000 | 1.000 [0.610, 1.000] | 1.000 |
| species | 6 | 0.333 | 1.000 [0.342, 1.000] | 0.333 |

## By `n50_bin`


**n50_bin = <10 kb** (n = 345)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 345 | 0.997 | 0.988 [0.970, 0.995] | 0.986 |
| class | 336 | 1.000 | 0.982 [0.962, 0.992] | 0.982 |
| order | 334 | 0.997 | 0.985 [0.965, 0.994] | 0.982 |
| family | 335 | 0.988 | 0.988 [0.969, 0.995] | 0.976 |
| genus | 342 | 0.985 | 0.988 [0.970, 0.995] | 0.974 |
| species | 345 | 0.594 | 0.971 [0.938, 0.987] | 0.577 |
