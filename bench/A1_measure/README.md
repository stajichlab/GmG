# A1 measurement (HPCC): peak memory and time of `popsketch run --from-sketches`

Checks the arithmetic in `../A1_memory_arithmetic/README.md`. Run: `REPO=$PWD sbatch bench/A1_measure/measure.sbatch` from a checkout on HPCC. Builds in `$SCRATCH`, makes symlink sets in `$SCRATCH` (deleted when the job ends), writes only `results.tsv`.

Sets: `n264` (C. neoformans), `n2002` (all seven SRA sets), `n4004` and `n6006` (the 2,002 sketches repeated 2x and 3x under new ids). `--flagged mark` keeps every sample in the analysis so n is the true matrix size; `--no-pca` versus default shows what the ordinations add. Replicas are exact copies (distance 0 to their original), so time at large n is not representative of structureless real data; memory is.

Expectation to falsify (from the arithmetic): peak RSS within 30% of 7 x 8 n^2 bytes plus the sketches. `results.tsv` columns: set, variant, n_input, in_analysis, max_rss_gb, wall_s, exit.
