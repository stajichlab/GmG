# Plan: a test framework for "popsketch at scale" (draft v3, reviewed twice)

## The question

Can we, for thousands of strains across several species: download each run, subsample to about 20-30x, sketch it at a good size/sensitivity trade-off, and give a browser (WASM) page the means to make tables, PCoA/PCA, trees, hypothesis tests and a query over the sketches? This is not the build plan (see `plan_web_sketches.md`). It is the set of measurements that tells us (1) what is possible and how well, and (2) what resources it takes, before we commit. v2 folds in an independent review (Fable) and a second (Opus) that read PopSeq; changes are marked *(review)* and *(review 2)*.

## Facts corrected or added

- The fetch script (`03_fetch_sketch.sbatch`) keeps every Nth read pair of the whole run, so cost scales with the run size. `step = ceil(bases/1e9)`, so it keeps between 0.5 and 1 Gb (a 1.01 Gb run keeps about 0.5 Gb): depth already varies twofold across the 2,002 sketches. `08_depth_series.sbatch` uses a fractional evenly spaced sampler instead. *(review)*
- "20-30x" must say which coverage. At 1 Gb, *Z. tritici* has about 18x k-mer coverage and about 28-55x in bases for 18-35 Mb genomes (my earlier "20-100x" was loose). The cap is per species (genome size x target), and the bases actually kept must be logged. *(review)*
- Sizes: a `.gsk` stores u64 hash + u32 abundance, 12 B per hash, and a read sketch keeps error k-mers below the cutoff (about 3 MB, not 1.2 MB). "5 B per hash" holds only for delta-varint hashes above the cutoff with abundances dropped (the compact format). *(review)*
- Code for distances, ordination and tree is in `crates/popsketch`, not `gmg-core`. `DistMat` holds dense n x n f64 for both `ani` and `mash`. Tree = neighbour-joining on mash (Jaccard); PCoA = containment ANI over the smaller set, which hides extra content. *(review)*
- Measured: 8-13 min wall on 4 cores per run at about 1 Gb (including download); native analysis about 1 min at 1,500 samples and 4.5 min at 3,000 (synthetic hashes).

## Working rules (owner, 2026-11-07)

1. **HPCC only.** All analysis and data processing runs on HPCC; the laptop is the orchestrator. The one exception (confirmed) is the browser timings in T5 and T7b, which need a real browser. `~/.s3cfg` is on both machines, so `s3cmd` (T8 and the upload) runs on HPCC too, from the owner's home; the credential is not copied anywhere else.
2. **Purgeable folders.** One-off and exploratory work lives under `~/projects/GmG/` on HPCC in named folders: `bench/<test-id>_<what>/` (T-numbered, safe to delete after the ledger is written), `scratch_notes/`, and `results/` only for what a document cites. Anything not under `bench/`, `analysis/` or `results/` is not kept.
3. **Temp files in `$SCRATCH`.** Downloads (.sra, FASTQ), subsamples and anything in-silico mixed are written to `$SCRATCH` (`/scratch/jstajich/`) inside the task and removed on exit; only sketches, traces and tables are kept. Nextflow `work/` holds only the declared small outputs.
4. **Code and tests.** Scripts under `analysis/` or `bench/`, Rust in `crates/popsketch` (and `popsketch-wasm`). Finished code is documented in `docs/popsketch.md`; each test above gets an automated check (a small fixture run in CI, a scale test marked `--ignored`, and the benchmark itself), not only a one-off run.
5. **PopSeq owns metadata.** Runs, accessions, genome sizes, labels and panels come from the PopSeq catalog (`/bigdata/stajichlab/shared/projects/Fungi/PopSeq`: `catalog/index.json`, `catalog/runs`, `catalog/genome_sizes.csv` with 113 species, `uv run popseq panel --min-depth ...`). GmG only reads them. If a field is missing or wrong, or another NCBI/SRA query is needed, the fix is proposed in PopSeq. Consequence: for the 2,002 existing sketches the manifest joins the existing accessions to PopSeq (they came from MATPredict tables, and a seeded `popseq panel` would pick different runs); `popseq panel` is for new species only. `01_select_runs.py` and `02_ena_urls.py` are retired once that join exists. The per-species base cap (genome size x target) reads `genome_sizes.csv` and `est_depth_x`; note PopSeq has *Z. tritici* at 39.7 Mb and *P. oryzae* at 41.0 Mb, so earlier 18-35 Mb figures for those do not apply.

**PopSeq issues** *(review 2)*: opened #23 (country and host) and #24 (library_strategy); the downloads log is the existing #14 (D3), and T1 is benchmark #12 (D1). Details: normalised `country` and `host` (free text: 316 distinct values for *Z. tritici*, including "", "missing" and two spellings of Eschikon; only 59% of runs have a country), `library_strategy` is WGS on every row so that filter does nothing, and the downloads log (T1). T9 and the web catalog must not re-parse these.

## Method

Each test has a question, a dataset, a metric, a pass criterion set **before** running, and a resource record. Order is cheap to expensive; a failed gate stops dependent tests. A failed gate changes the plan; the criterion is not relaxed afterwards.

**Execution layer: Nextflow** on HPCC (modules exist). The head job runs on `epyc` or `batch` (30-day limit), not `short`. One task per run reads a manifest and runs under `scratch '$SCRATCH'` (single quotes; `$SCRATCH` is set only inside a job and is node-local NVMe of about 1.5 TB shared by every job on that node, so the budget is per node, not 1 TB overall). Without the `scratch` directive `work/` keeps everything a task writes, so reads must stay in `$SCRATCH`; `work/` then holds sketches and logs (about 40 GB at 8 MB x 5,000, about 5 GB compact) under `~/projects`, not the shared folder (`/bigdata` is 94% full). `short` is 2 h: large runs use `errorStrategy 'retry'` that moves the task to `epyc` with a longer time limit. `-resume` gives "what has run / what is left"; `maxForks` and `queueSize` keep NCBI and the queue polite; failures follow PopSeq's states: a 404 at NCBI is `unavailable_at_ncbi` (not transient, never retried at NCBI) and the run is routed to ENA if the catalog says it is there; the route is chosen up front from the catalog's `ncbi_status` (`in_ncbi` 220,296; `ena_only` 621; `not_checked` 24,905). Only runs that fail every route go to `DATA_PROBLEMS.md` and `analysis/popgen_sra/excluded_runs.tsv`. Laptop needs under 1 GB. **The resource record is Nextflow's `trace.txt`** (wall, CPU, peak RSS, bytes read/written, exit status), plus browser timings for T5/T7b; costs are fitted (cost = a + b x run size), not single points. Bases actually kept and k-mer coverage are logged per sample.

**Data.** Existing: the 7 SRA sets (2,002 sketches), the Coccidioides and *B. dendrobatidis* sets, the 12-run deep sets (`08_select_deep.py`). Needed: a larger-genome species (100-300 Mb). Synthetic data only for pure scaling in n.

**Truth.** *(review)* No test judges itself by a deep run of the same pipeline: the old 6 Gb reference had a failed peak estimate for 9 of 12 *C. neoformans* runs (before #157) and is shallower than 30x for a 300 Mb genome. Use outside truth where it exists: SNP trees and distances for Cocci (446-tip tree) and Bd (VCF), public VCFs for *Z. tritici* and *C. neoformans*, scored with `--compare-tree` and `--compare-pcs`. Freeze one noise floor, measured at the reference depth, so shallow depths cannot pass by having a bigger floor. *(review 2)* **The truth must cover the T2 samples.** `08_select_deep.py` picks the 12 smallest 7-20 Gb runs per panel, which are in no VCF, and 12 runs give only 66 pairs (too few for ARI or NN agreement). So T2 runs are chosen by joining PopSeq accessions to the truth sets first: Cocci and Bd are the primary arms; *Z. tritici* and *C. neoformans* only for strains that are in a public VCF. Minimum n (stated before running): 30 runs per arm, 435 pairs.

## Order (most learned per core-hour) *(review)*

1. **Analytic and no-HPCC:** A1 (memory and time arithmetic), A2 (sampling variance and saturation), T8 (hosting throwaway file).
2. **T2 + T10 + T11** on the already listed deep runs, with the SNP truth and in-silico mixtures.
3. **T4, T5, T6, T9, T7b** on the results.
4. **T1**, cut down, last on HPCC. **T7** after the T6 decision.

## Tests

### A1. Memory and time arithmetic (before any benchmark)
**Measured 2026-11-08 (`bench/A1_memory_arithmetic/README.md`, `bench/A1_measure/`):** the matrix arithmetic below is right but small; peak RSS was 17 GB at 2,002 samples (about 7-9 MB per sample, the sketches themselves), distances dominate time, and the apparent cliff above about 2,000 pooled samples was an artefact of the replicated test sets (D1b: the real pooled set of 1,822 samples used the bitset with 2.8 M variable hashes, 0.59 GiB). Implications: precomputed per-species matrices for the page, group by species, drop unused split halves.
`DistMat` is 2 n² f64, plus a subset copy, `one_minus_ani`, and squared and centred PCoA matrices: about 7 n² f64 arrays, roughly 5.6 GB at 10k samples (not 800 MB). NJ is cubic: 1.3 s at 1,500 and 13 s at 3,000 projects to about 8 min at 10k. Compute n_max for a 4 GB / 60 s budget on paper, then refactor candidates (condensed f32 triangle, one distance only, RapidNJ or a faster NJ, partial/randomized eigen) and benchmark those.

### D1. Keep n small enough for the bitset path: collapse near-identical samples *(owner decision 2026-11-08)*
Design: analyse **per species**; before ordination and tree, collapse samples whose ANI to a representative is at or above a threshold, keep one representative per cluster and carry the cluster size as a weight, so nothing is silently dropped. The cap on n is then set by the bitset limit (`n x V / 8 <= BITSET_MAX_BYTES`, now 4 GiB) on the representatives, not by raw sample count. Also report the distance method and V in the `distances_all` stage line, and make the limit a flag. Label tests must use the weights (clonal samples are not independent replicates).

**Measured (`bench/D1_collapse/`, HPCC job 31876237, 95 s, existing sketches, flagged samples excluded so n is below 300):** representatives left by leader clustering (every member within the threshold of its representative); `components` is the chained upper bound on collapse.

| set (n) | ANI 0.99 | 0.995 | 0.999 | 0.9995 | 0.9999 |
|---|---|---|---|---|---|
| z_tritici (259) | 160 | 187 | 197 | 203 | 236 |
| c_neoformans (268) | 5 | 21 | 61 | 71 | 175 |
| c_beticola (282) | 1 | 8 | 168 | 173 | 202 |
| p_nodorum (265) | 1 | 86 | 196 | 199 | 231 |
| b_cinerea (294) | 1 | 34 | 179 | 187 | 250 |
| n_discreta (203) | 9 | 11 | 119 | 166 | 203 |
| p_oryzae (251) | 6 | 12 | 124 | 182 | 245 |

What it says: only the clonal sets collapse a lot at a tight threshold (C. neoformans 268 to 61 at 0.999; P. oryzae and N. discreta about half). The other sets keep 60-80% of samples at 0.999, so collapse alone will not make a 10,000-strain species fit; the n cap still has to bite. Thresholds of 0.995 and below merge whole populations (c_beticola, p_nodorum, b_cinerea go to 1-34 clusters), so they are not usable for collapse. **Not yet answered:** what threshold is distinguishable from sketch noise at the capped depth. At 0.9999 most samples stay separate, which suggests the noise floor, not biology, sets the usable threshold. Next measurement: take the split-half ANI of each sample (already in the sketches) as the floor and choose the threshold from it; then recount, and measure pooled V for the representatives to turn the cap into a number.

**D1b (`bench/D1_collapse/results_floor.tsv`, `pooled.tsv`, job 31880606): thresholds from the noise floor, and pooled runs.** The split-half floor (90th percentile of 1 - ANI) is 0.0025-0.0046 per species, and is pessimistic (each half has half the coverage). Representatives left at 1 floor: z_tritici 187 of 259, c_neoformans 41 of 268, c_beticola 21 of 282, p_nodorum 161 of 265, b_cinerea 143 of 294, n_discreta 12 of 203, p_oryzae 47 of 251. At 3 floors or more nearly everything merges, so 1-2 floors is the usable range, and one floor already merges most of the structure of the low-diversity sets (c_beticola, n_discreta). Because the floor is pessimistic, the threshold should be chosen per species with the floor in view, not fixed. Pooled runs (`--no-pca`, all on the bitset path, limit 4 GiB):

| pooled set | samples | variable hashes | presence matrix | wall |
|---|---|---|---|---|
| everything (1,554 in analysis) | 1,822 | 2.79 M | 0.59 GiB | 59 s |
| representatives at 1 floor | 612 | 2.44 M | 0.17 GiB | 18 s |
| at 2 floors | 235 | 1.79 M | 0.05 GiB | 5 s |
| at 3 floors | 104 | 1.31 M | 0.02 GiB | 3 s |

What it changes: the bitset limit is far from these set sizes, so the n cap can be set by what the plots can show, not by the 4 GiB limit; the limit is about n = 12,000 at the current V (an upper estimate, since V grows). Collapse is still useful for clarity and for the label tests. Note `--max-samples` defaults to 1,500 and refuses larger sets. Code added: `dist::compute_planned` reports the method and V (printed on the `distances_all` stage line), and `--bitset-max-gib` sets the limit.

**D1c (`bench/D1_collapse/species_threshold.tsv`, job 31884111): the per-species threshold.** The floor from 20 samples was replaced by the split-half distance of every sample (the halves are on disk for all of them), and each sample's distance to its nearest neighbour was added.

| species (n) | median split-half | p90 split-half | nearest-neighbour distance p10 / p50 / p90 | reps at 0.5 x median | at median | at p90 |
|---|---|---|---|---|---|---|
| z_tritici (259) | 0.0019 | 0.0048 | 0.0001 / 0.0116 / 0.0136 | 198 | 192 | 188 |
| c_neoformans (268) | 0.0003 | 0.0017 | 0.00003 / 0.0001 / 0.0033 | 149 | 94 | 51 |
| c_beticola (282) | 0.0014 | 0.0034 | 0.00004 / 0.0002 / 0.0028 | 171 | 156 | 63 |
| p_nodorum (265) | 0.0022 | 0.0046 | 0.0001 / 0.0027 / 0.0055 | 196 | 188 | 124 |
| b_cinerea (294) | 0.0014 | 0.0029 | 0.0001 / 0.0003 / 0.0032 | 182 | 177 | 133 |
| n_discreta (203) | 0.0022 | 0.0039 | 0.0002 / 0.0009 / 0.0022 | 128 | 59 | 12 |
| p_oryzae (251) | 0.0022 | 0.0040 | 0.0001 / 0.0006 / 0.0018 | 118 | 66 | 19 |

Reading: the species fall in two kinds. In z_tritici and p_nodorum the typical sample has no relative within noise (nearest-neighbour p50 is 0.003-0.012, above the floor), so collapse removes only 20-30% whatever the rule; these are recombining populations and collapse is not the tool. In c_neoformans, c_beticola, b_cinerea, n_discreta and p_oryzae the typical sample has a relative at or below the floor (clonal or outbreak sets), and the rule decides the result. The p90 rule collapses n_discreta to 12 and p_oryzae to 19, which discards structure that is probably real (the p90 is set by a few low-coverage samples).

**Rule (proposed, not yet validated): collapse threshold = 1 - median split-half distance of that species, computed from all samples that carry halves, with the p90 shown as the looser alternative.** Pairs closer than the typical technical-replicate distance cannot be told apart by this sketch, so one representative stands for them. The halves have half the coverage, so the real noise between full sketches is lower and the rule errs toward merging a little too much. **Validation still to do:** compare a full sketch with its own lower-depth subsamples (the `depth/` sets hold the same runs at several depths) to see how the real distance of a sample to itself depends on coverage; a per-pair rule using each sample's own split-half distance is the refinement if the median proves too blunt.

**D1d (`bench/D1_collapse/depth_floor.tsv`, `depth_self.tsv`, job 31890005): depth check.** 12 runs each of z_tritici and c_neoformans, subsampled by stride to 1, 2 and 4 Gb (about 25, 50 and 100x for z_tritici; 52, 105 and 210x for c_neoformans). Median split-half distance by depth: z_tritici 0.00065 / 0.00016 / 0.00006; c_neoformans 0.00033 / 0.00011 / 0.00004. So the noise falls about 3-4 fold per doubling of depth, not by the square root. A run's distance to its own lower-depth subsample is at most 0.0001 in nearly every case (medians 0 to 0.00006; one c_neoformans run reached 0.001 between 1 and 2 Gb).

What it does and does not show. It does not validate the rule: the 1, 2 and 4 Gb subsamples are nested draws from the same reads, and containment between nested draws is biased toward agreement, so those distances are a lower bound on noise, not independent replicates (the split-half distance, built from alternate reads, is the honest one). It does show that the floor is mostly a coverage effect: the 12 high-depth runs have a median floor about 3 times lower at 1 Gb (0.00065) than the whole z_tritici set (0.0019), which includes runs with fewer than 1 Gb. A per-species median therefore mixes coverage differences with biology. **Revised proposal: a per-pair rule, with the pair threshold taken from the two samples' own split-half distances (for example the mean of the two), reported with the per-species median as the summary.** That needs the per-sample split-half distance written out (not yet in the outputs), so the next step is to expose it and recount. Doubling the cap from 1 to 2 Gb would cut the floor about 4 fold at twice the sketch size; that trade belongs with T2. Caveat: 12 runs per species, two species.

**D1e (`bench/D1_collapse/per_pair.tsv`, `split_half_<species>.tsv`, job 31895509): the per-pair rule.** `popsketch run` now writes `split_half.tsv` (each sample's own split-half distance). Representatives left (leader clustering in file order; a sample joins a representative when 1 - ANI is at most the threshold):

| species (n) | species median (D1c) | per-pair mean | per-pair max |
|---|---|---|---|
| z_tritici (259) | 192 | 191 | 190 |
| c_neoformans (268) | 94 | 86 | 77 |
| c_beticola (282) | 156 | 146 | 138 |
| p_nodorum (265) | 188 | 184 | 176 |
| b_cinerea (294) | 177 | 167 | 160 |
| n_discreta (203) | 59 | 52 | 49 |
| p_oryzae (251) | 66 | 54 | 39 |

The species-median column reproduces D1c exactly, so the runs agree. The per-pair mean removes 0-18% more representatives than the species median (largest in p_oryzae, 66 to 54, and n_discreta, 59 to 52); the max rule goes further (p_oryzae 39). Largest clusters grow with the per-pair rules (c_neoformans 34 to 93, p_oryzae 55 to 123). So the choice changes the count modestly and the cluster sizes more; it does not rescue or break any set. **Decision proposed: per-pair mean** (each pair's threshold is the mean of the two samples' own split-half distances; samples without halves use the species median), because it follows each sample's coverage and needs no species-wide constant, with the species median kept as the reported summary. Still unvalidated against independent replicates: no set here has real technical replicates of the same sample, which is the test that would say whether the noise estimate is right. Leader clustering depends on input order; the count is one valid clustering, not the minimum. **Next in code:** a collapse step in `popsketch` using this rule, writing representatives and cluster sizes, and the label tests using the sizes as weights.

### A2. FracMinHash validity (predict first)
At scaled 200 a 20 Mb genome gives about 100k hashes; sampling SD of 1 - ANI near 0.001 is about 2e-5, far below the measured floors (0.001-0.0025). So floors come from reads, errors and the cutoff, not from sampling. At the other end, below about 80% ANI, C^(1/21) rests on a handful of shared hashes (saturation). Predict both regimes; run only scaled 200 and 1000 at k=21 unless the prediction fails (so k=17/25 are dropped).

### T2. Subsample design and depth
From one download per run: fractional even sampler at 5, 10, 15, 20, 30, 50x (bases and k-mer coverage both stated); PopSeq's default `windows` subsample (N spots from 8 evenly spaced windows, `fastq-dump -N/-X`; spec 002 forbids `head` for cross-sample measurements, so first-N is a negative control only); a random subset. Arms: Cocci, Bd, and the *Z. tritici*/*C. neoformans* strains with VCFs, at least 30 runs each. The large-genome arm: PopSeq has no species near 300 Mb. Over 100 Mb it has only *Blumeria graminis* (141 Mb, 959 runs), *Puccinia triticina* (123, 610), *Ophiocordyceps sinensis* (111, 238), *Melampsora laricis-populina* (101, 267); these are obligate parasites with host reads, so that arm is treated as a T10 case as well as a depth case, and the cost figures use about 140 Mb x 25 = 3.5 Gb, not 9 Gb. Species is chosen when T2 is scheduled (default *B. graminis*). Measure whether `windows` really transfers less than the whole run (GmG's reply to 004, item 6, says it does not); the reads are kept for the job (the 03 script deletes them). Metrics against outside truth and the frozen floor: Spearman correlation of pairs above the floor, mean difference versus the frozen floor, ARI of clusters, and nearest-neighbour agreement **counted only where the gap between the first and second neighbour exceeds the floor** (near-clonal sets like Bd, where 54% of pairs are below the floor, and *Z. tritici* at 83-92% between 1 and 6 Gb, cannot reach a flat 95%). Output: per-species cap rule and the T1 cost saved.

### T3. Compact size (folded into A2)
A2 predicts the k/scaled behaviour. T3 keeps only the measurement: bytes per sample in the compact format at scaled 200 and 1000 (assumptions stated above), and clade recovery (Cocci, Bd). Other k and scaled only if A2's prediction fails.

### T4. Native scaling in n
After A1: real sketches resampled to n = 500, 1k, 2k, 5k, 10k, for the current code and the refactored candidates. Time and RSS per stage. Pass: stated n_max for 4 GB / 60 s.

### T5. WASM and the page
Chrome, Firefox, Safari, desktop and one phone, in a worker. Load and parse time, compute time, peak memory, native/wasm ratio, bytes fetched; real OOM n (wasm32 ceiling). Parity with native at every n. **Design decision made now:** the page fetches precomputed condensed matrices and ordinations (including the hash-PCA coordinates, which cannot be rebuilt from a matrix and are the default ordination) per species (about 50 MB at 5k in f32); sketches (about 3 MB per sample in a bundle, about 15 GB at 5,000) are loaded only for a query or added samples. Pass: time to first plot.

### T6. Cross-species behaviour
Pool the 7 species plus distant controls: do 1 - ANI and containment separate within-species, between-species and between-genus, and where do they saturate. Pass: a written rule for when two samples are comparable (expected: group by species before ordinating).

### T7. Query tier
Hash-to-sample inverted index for 2,002 and replicated 5,000 samples: size, build time, "closest samples to this sketch" and "contains this sequence" latency versus brute force; sharded by hash prefix (no Range on the host). Pass: under 2 s at 5,000 with a per-species index small enough to hold. Waits on T6.

### T7b. Browser-side database
Whole-file `fetch` from S3 (CORS only, no Range) into DuckDB-wasm (`registerFileBuffer`, Parquet catalog) and into Rust/WASM (sorted delta-coded inverted index shard), thrown away on tab close. 300 and 1,000 samples per species: memory, load time, query latency, desktop and phone. About 90k hashes per *C. neoformans* sample means 27M rows at 300 samples; cross-species tables are too large, so shard per species or hash prefix. Expectation, not measured: DuckDB for metadata, Rust for the hash index.

### T8. Hosting and transfer
NRP S3 `s3://stajichlab/sketch/fungi/` (owner's `~/.s3cfg`, on laptop and HPCC, `s3cmd -P`), throwaway file first: CORS from `https://gmg.fungalgenomes.org`, exposed `ETag`/`Content-Length`, Range support, concurrent GET throughput, egress limits. Pass: first view of a precomputed species under 30 s on a good connection (met by T5's design, not by fetching sketches).

**T8 partial results (2026-11-09, from the laptop; NRP S3 came back after maintenance).** A throwaway folder `s3://stajichlab/sketch/fungi/_t8_throwaway/` (8 MiB random file and a small JSON, uploaded with `s3cmd -P`) is served anonymously at `https://s3-west.nrp-nautilus.io/stajichlab/sketch/fungi/_t8_throwaway/<file>` (path-style URL; HTTP/2). Passes: anonymous GET and HEAD work over HTTPS; `Accept-Ranges: bytes` and a Range request returns `206` with the right `Content-Range`; `ETag` and `Content-Length` are returned. **Fails: CORS.** The bucket has no CORS configuration (`s3cmd info`), a GET with `Origin: https://gmg.fungalgenomes.org` returns no `Access-Control-Allow-Origin`, and the preflight `OPTIONS` returns 403, so the page on `gmg.fungalgenomes.org` cannot read these objects until a CORS rule is set on the bucket (S3 CORS is per bucket, not per prefix). Proposed rule, needs the owner's go-ahead: origins `https://gmg.fungalgenomes.org` (and `http://localhost:8000` for development), methods `GET` and `HEAD`, allowed headers `Range` (and `If-None-Match`), exposed headers `ETag`, `Content-Length`, `Content-Range`, `Accept-Ranges`, max age 3600; undo with `s3cmd delcors`. Throughput from the laptop (one location, one afternoon, noisy): one GET 3.4 MB/s; 4 parallel 2.4 MB/s aggregate; 8 parallel 6.4 MB/s aggregate, so parallel GETs do not scale per connection from here; a 3 GB species bundle would take 8-20 minutes and must not be the first view, whereas a precomputed matrix of 300 samples (about 0.4 MB as float32) is immediate. **CORS rule set (owner's go-ahead, 2026-11-09):** the rule above was applied to `s3://stajichlab` with `s3cmd setcors` (origins `https://gmg.fungalgenomes.org` and `http://localhost:8000`, GET and HEAD, allowed headers `Range` and `If-None-Match`, exposed `ETag`, `Content-Length`, `Content-Range`, `Accept-Ranges`, max age 3600; undo with `s3cmd delcors`). Checked with `curl`: a GET and a Range GET with the allowed origin return `Access-Control-Allow-Origin` for that origin and the exposed headers; the preflight (`Access-Control-Request-Headers: range,if-none-match`) returns 200 with the allowed headers; `http://localhost:8000` is allowed; another origin (`https://example.org`) gets no `Access-Control-Allow-Origin`; a request without `Origin` is unchanged. **Not yet done: the same read from a real browser tab on `gmg.fungalgenomes.org`** (the browser extension was not connected), and egress limits, browser throughput and throughput from another network. The throwaway files were deleted after the checks; recreate them with the T8 commands to repeat.

### T9. Hypothesis tests
*(review)* No PERMANOVA exists in the code yet. Build it with PERMDISP (unequal spread is mistaken for a location effect when groups are unbalanced, as in *P. nodorum*). Permute within study (studies are clustered and often clonal, so samples are not exchangeable), not all labels; add planted-effect sets to measure power. Pass: calibrated p-values and a confounding warning.

### T10. Heterozygosity, ploidy and contamination at the capped depth *(review, new)*
At 20-30x, heterozygous k-mers sit near 10-15x, close to the abundance cutoff, and the shoulder test needs about 12x k-mer coverage. Mix reads in silico: a contaminant at 1/3/10/30% (a second species), and a diploid such as *C. albicans*, at each depth. Report both distances (the tree uses mash, PCoA uses containment ANI, so they disagree under extra content) and flag recall (`size_high`, `abundance_shoulder`, `unrelated_to_cohort`). May change go/no-go for the depth cap.

### T11. Reads versus assemblies *(review, new)*
The database will mix SRA read sketches and assemblies. Take strains with both (BFD/Cocci sets) and measure the offset in both distances; decide whether they may share an ordination. Pass: the read-versus-assembly offset is below the frozen floor.

### T1. Acquisition cost (cut down, last) *(review 2)*
This is PopSeq benchmark 7.1, not a new design. Route 0 in PopSeq spec 002 is AWS Open Data over HTTPS (`curl -C -`, median 10.3 MB/s over 4,211 runs), then `prefetch`, then ENA, with `fasterq-dump` as the last resort. Run 20 runs stratified by size through the PopSeq routes plus the current `parallel-fastq-dump`, timed at several hours of the day, and write rows in PopSeq's log format (`benchmarks/downloads.csv`, which does not exist yet; contribute it there). The Nextflow task adds a mandatory mate-count check, but **not** against ENA `read_count`: for 172 of 200 sampled ERR runs ENA's `read_count` counts reads, not pairs (SRR 0 of 200, DRR 0 of 200; `benchmarks/2026-11-07_read_count_convention.csv` in PopSeq). The rule PopSeq adopted (spec 002, PR #25): verify SRA-derived FASTQ against NCBI `spots_with_mates` (each mate must equal it), and ENA-sourced FASTQ against ENA's MD5 plus equal mate counts. `03_fetch_sketch.sbatch` has no check, so a truncated dump would silently give a shallower sketch. Unpaired runs: 27 of 200 sampled ERR runs have `spots_with_mates = 0` at NCBI, so `parallel-fastq-dump --split-files` cannot produce mates for them; ENA has proper `_1`/`_2` files, so the route for such runs is ENA (decided from the catalog before fetching), and the task records the layout it actually got. Until PopSeq #13 lands, the check is done in the task with `spots_with_mates` from the run's NCBI metadata; afterwards it calls PopSeq's verifier. Measured on compute nodes so far (PopSeq, n = 9-12 per route, early): AWS 10.3 MB/s, `prefetch` 9.5 MB/s, ENA HTTPS 0.88 MB/s; Aspera does not work (no key; `prefetch -t fasp` fails). ENA truncation is dropped (PopSeq settled it). Record wall, CPU, bytes in, scratch peak (a `du` sampler inside the task; `trace.txt` has no scratch column), failure state. A 140 Mb genome at 25x keeps about 3.5 Gb.

## Resource ledger (to fill in)

| item | per sample | per 5,000 samples | source |
|---|---|---|---|
| download+subsample core-hours (0.5-0.9 at 1 Gb) | T1/T2 | | trace.txt |
| scratch peak | T1 | | trace.txt |
| compact sketch bytes | T3 | | |
| analysis time and RAM, native / wasm | A1, T4 / T5 | | |
| browser memory | T5, T7b | | |
| storage, egress | T8 | | |
| index | T7 | | |

Budgets are separate: T2/T10 about 36 deep runs of reads kept (core-hour estimate to be made from T1's fit, not assumed); T1 20 runs x 2 methods; T4 small; T5-T9 on the laptop. Estimates only until A1 and the first trace exist. Under 1 TB transient scratch; run when HPCC load is low.

## Decisions the results feed

Depth cap rule and subsample method; k, scaled, file format; per-species n_max (native and browser) and whether a server tier is needed beyond it; static sharded index or service; reads-and-assemblies policy; host and cost per 1,000 views.

## Risks

Truth only exists for some species; small numbers of runs per species; one large genome may not represent others; NCBI/ENA throughput varies by hour; browser memory differs by device; T6 may show the sketch cannot compare beyond genus, which limits the multi-species view to a catalog plus search.
