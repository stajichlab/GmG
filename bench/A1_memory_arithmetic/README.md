# A1: memory and time arithmetic for `popsketch run` (paper estimate, 2026-11-08)

Purpose: work out n_max for a memory and time budget before benchmarking. Purgeable: nothing here is data; it is arithmetic plus a pointer to what to measure on HPCC (rule 1: measurements run there).

## Dense n x n f64 arrays alive during a run (`crates/popsketch/src`)

| array | where |
|---|---|
| `ani` and `mash` of all candidates (`d_all`) | `dist.rs:79-80` (`DistMat`) |
| `ani` and `mash` of the included subset (`d_all.subset(...)`) | `pipeline.rs:779` |
| `d1 = one_minus_ani()` | `pipeline.rs:787`, `dist.rs:22` |
| squared distances `d2` | `ordination.rs:324` (PCoA) |
| double-centred `b` | `ordination.rs:329` (PCoA) |
| Gram matrix `g` | `ordination.rs:300` (hash PCA, `gram_matrix`) |
| UPGMA and NJ working copies | `cluster.rs`, `tree.rs:23` (to check) |

A run can hold about 7 of these at once (more with UPGMA/NJ copies): `d_all` (2) + subset (2) + `d1` (1) + `d2` and `b` (2). `top_eigs` (block iteration, `ordination.rs:130`) avoids Jacobi's extra n x n, but Jacobi (`ordination.rs:26`) allocates one more.

## Memory (one f64 array = 8 n^2 bytes)

| n | one array | 7 arrays | 9 arrays | one f32 condensed (2 n^2) |
|---|---|---|---|---|
| 1,500 | 18 MB | 126 MB | 162 MB | 4.5 MB |
| 3,000 | 72 MB | 504 MB | 648 MB | 18 MB |
| 5,000 | 200 MB | 1.4 GB | 1.8 GB | 50 MB |
| 10,000 | 800 MB | 5.6 GB | 7.2 GB | 200 MB |

n_max for 7 live f64 arrays (56 n^2 bytes): 2 GB gives n of about 6,000; 4 GB gives about 8,400. The wasm32 address space is 4 GB, and the browser also holds the sketch data, so the practical browser limit is lower (T5 measures it).

## Time

Neighbour-joining is cubic. Measured on structureless synthetic hash sets (`docs/popsketch.md`): 1.3 s at 1,500, 13 s at 3,000 (the measured exponent is about 3.3). Scaling the 3,000 point by (n/3000)^3: about 60 s at 5,000 and about 480 s (8 min) at 10,000. Hash PCA and PCoA were about 125 s each at 3,000 on data with no structure, the eigen-solver's slowest case; their growth in n is not extrapolated here (T4 measures it).

## What to change, in order of payoff

1. One distance only per consumer: tree on mash needs `mash`; clusters and PCoA need `1 - ani`. Do not keep both dense after the subset.
2. Condensed f32 upper triangle (2 n^2 bytes instead of 8 n^2 per array). Check that f32 keeps differences smaller than the noise floor (0.001-0.0025) and the tie-breaking in UPGMA: f32 has about 7 digits, so distances near 0.001 keep 1e-7 resolution, fine for these floors but this must be a test, not an assumption.
3. Drop the subset copy when nothing was excluded.
4. A faster NJ (RapidNJ or canonical-search variants) if 10k trees are wanted in the browser; otherwise tree only for n of about 3,000 or fewer.

## To measure on HPCC (T4)

`/usr/bin/time -v` around `popsketch run --from-sketches` on one real set (C. neoformans 300 samples, and replicated sets at 1k, 2k, 5k) with `--no-pca` on and off; compare peak RSS with the table above. Expectation to be falsified: peak RSS is within 30% of 7 x 8 n^2 plus the sketches.
