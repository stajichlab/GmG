# Plan: a test framework for "popsketch at scale" (draft v2, reviewed)

## The question

Can we, for thousands of strains across several species: download each run, subsample to about 20-30x, sketch it at a good size/sensitivity trade-off, and give a browser (WASM) page the means to make tables, PCoA/PCA, trees, hypothesis tests and a query over the sketches? This is not the build plan (see `plan_web_sketches.md`). It is the set of measurements that tells us (1) what is possible and how well, and (2) what resources it takes, before we commit. v2 folds in an independent review (Fable) that checked this plan against the code; changes are marked *(review)*.

## Facts corrected or added

- The fetch script (`03_fetch_sketch.sbatch`) keeps every Nth read pair of the whole run, so cost scales with the run size. `step = ceil(bases/1e9)`, so it keeps between 0.5 and 1 Gb (a 1.01 Gb run keeps about 0.5 Gb): depth already varies twofold across the 2,002 sketches. `08_depth_series.sbatch` uses a fractional evenly spaced sampler instead. *(review)*
- "20-30x" must say which coverage. At 1 Gb, *Z. tritici* has about 18x k-mer coverage and about 28-55x in bases for 18-35 Mb genomes (my earlier "20-100x" was loose). The cap is per species (genome size x target), and the bases actually kept must be logged. *(review)*
- Sizes: a `.gsk` stores u64 hash + u32 abundance, 12 B per hash, and a read sketch keeps error k-mers below the cutoff (about 3 MB, not 1.2 MB). "5 B per hash" holds only for delta-varint hashes above the cutoff with abundances dropped (the compact format). *(review)*
- Code for distances, ordination and tree is in `crates/popsketch`, not `gmg-core`. `DistMat` holds dense n x n f64 for both `ani` and `mash`. Tree = neighbour-joining on mash (Jaccard); PCoA = containment ANI over the smaller set, which hides extra content. *(review)*
- Measured: 8-13 min wall on 4 cores per run at about 1 Gb (including download); native analysis about 1 min at 1,500 samples and 4.5 min at 3,000 (synthetic hashes).

## Working rules (owner, 2026-11-07)

1. **HPCC only.** All analysis and data processing runs on HPCC; the laptop is the orchestrator. Exceptions that need a laptop (to be confirmed with the owner): browser timings in T5 and T7b, and the S3 upload and CORS check in T8 (credentials are on the laptop, never on HPCC).
2. **Purgeable folders.** One-off and exploratory work lives under `~/projects/GmG/` on HPCC in named folders: `bench/<test-id>_<what>/` (T-numbered, safe to delete after the ledger is written), `scratch_notes/`, and `results/` only for what a document cites. Anything not under `bench/`, `analysis/` or `results/` is not kept.
3. **Temp files in `$SCRATCH`.** Downloads (.sra, FASTQ), subsamples and anything in-silico mixed are written to `$SCRATCH` (`/scratch/jstajich/`) inside the task and removed on exit; only sketches, traces and tables are kept. Nextflow `work/` holds only the declared small outputs.
4. **Code and tests.** Scripts under `analysis/` or `bench/`, Rust in `crates/popsketch` (and `popsketch-wasm`). Finished code is documented in `docs/popsketch.md`; each test above gets an automated check (a small fixture run in CI, a scale test marked `--ignored`, and the benchmark itself), not only a one-off run.
5. **PopSeq owns metadata.** Runs, accessions, genome sizes, labels and panels come from the PopSeq catalog (`/bigdata/stajichlab/shared/projects/Fungi/PopSeq`: `catalog/index.json`, `catalog/runs`, `catalog/genome_sizes.csv` with 113 species, `uv run popseq panel --min-depth ...`). GmG only reads them. If a field is missing or wrong, or another NCBI/SRA query is needed, the fix is proposed in PopSeq. Consequence: `analysis/popgen_sra/01_select_runs.py` and `02_ena_urls.py` duplicate PopSeq and should be replaced by a manifest built from `popseq panel`; the per-species base cap (genome size x target) reads `genome_sizes.csv`.

## Method

Each test has a question, a dataset, a metric, a pass criterion set **before** running, and a resource record. Order is cheap to expensive; a failed gate stops dependent tests. A failed gate changes the plan; the criterion is not relaxed afterwards.

**Execution layer: Nextflow** on HPCC (modules exist; run the head job there). One task per run reads a manifest (`runs.tsv`), works in node `$SCRATCH` (about 1 TB transient is fine), and declares only the sketch as output, so `work/` holds sketches and logs (about 40 GB at 8 MB x 5,000, about 5 GB compact) under `~/projects`, not the shared folder. `-resume` gives "what has run / what is left"; `maxForks` and `queueSize` keep NCBI and the queue polite; permanent failures (404) go to a failure list that feeds `DATA_PROBLEMS.md` and `analysis/popgen_sra/excluded_runs.tsv`. Laptop needs under 1 GB. **The resource record is Nextflow's `trace.txt`** (wall, CPU, peak RSS, bytes read/written, exit status), plus browser timings for T5/T7b; costs are fitted (cost = a + b x run size), not single points. Bases actually kept and k-mer coverage are logged per sample.

**Data.** Existing: the 7 SRA sets (2,002 sketches), the Coccidioides and *B. dendrobatidis* sets, the 12-run deep sets (`08_select_deep.py`). Needed: a larger-genome species (100-300 Mb). Synthetic data only for pure scaling in n.

**Truth.** *(review)* No test judges itself by a deep run of the same pipeline: the old 6 Gb reference had a failed peak estimate for 9 of 12 *C. neoformans* runs (before #157) and is shallower than 30x for a 300 Mb genome. Use outside truth where it exists: SNP trees and distances for Cocci (446-tip tree) and Bd (VCF), public VCFs for *Z. tritici* and *C. neoformans*, scored with `--compare-tree` and `--compare-pcs`. Freeze one noise floor, measured at the reference depth, so shallow depths cannot pass by having a bigger floor.

## Order (most learned per core-hour) *(review)*

1. **Analytic and no-HPCC:** A1 (memory and time arithmetic), A2 (sampling variance and saturation), T8 (hosting throwaway file).
2. **T2 + T10 + T11** on the already listed deep runs, with the SNP truth and in-silico mixtures.
3. **T4, T5, T6, T9, T7b** on the results.
4. **T1**, cut down, last on HPCC. **T7** after the T6 decision.

## Tests

### A1. Memory and time arithmetic (before any benchmark)
`DistMat` is 2 n² f64, plus a subset copy, `one_minus_ani`, and squared and centred PCoA matrices: about 7 n² f64 arrays, roughly 5.6 GB at 10k samples (not 800 MB). NJ is cubic: 1.3 s at 1,500 and 13 s at 3,000 projects to about 8 min at 10k. Compute n_max for a 4 GB / 60 s budget on paper, then refactor candidates (condensed f32 triangle, one distance only, RapidNJ or a faster NJ, partial/randomized eigen) and benchmark those.

### A2. FracMinHash validity (predict first)
At scaled 200 a 20 Mb genome gives about 100k hashes; sampling SD of 1 - ANI near 0.001 is about 2e-5, far below the measured floors (0.001-0.0025). So floors come from reads, errors and the cutoff, not from sampling. At the other end, below about 80% ANI, C^(1/21) rests on a handful of shared hashes (saturation). Predict both regimes; run only scaled 200 and 1000 at k=21 unless the prediction fails (so k=17/25 are dropped).

### T2. Subsample design and depth
From one download per run: fractional even sampler at 5, 10, 15, 20, 30, 50x (bases and k-mer coverage both stated); first-N reads; a random subset. 12 runs each of *C. neoformans*, *Z. tritici* and one genome >= 100 Mb; the reads are kept for the job (the 03 script deletes them). Metrics against outside truth and the frozen floor: Spearman correlation of pairs above the floor, mean difference versus the frozen floor, ARI of clusters, and nearest-neighbour agreement **counted only where the gap between the first and second neighbour exceeds the floor** (near-clonal sets like Bd, where 54% of pairs are below the floor, and *Z. tritici* at 83-92% between 1 and 6 Gb, cannot reach a flat 95%). Output: per-species cap rule and the T1 cost saved.

### T3. Sketch parameters
Scaled 200 and 1000, k=21, same reads. Bytes per sample in the compact format (stated assumptions, see above), and resolution on known close pairs and clade recovery (Cocci, Bd). Other k and scaled only if A2's prediction fails.

### T4. Native scaling in n
After A1: real sketches resampled to n = 500, 1k, 2k, 5k, 10k, for the current code and the refactored candidates. Time and RSS per stage. Pass: stated n_max for 4 GB / 60 s.

### T5. WASM and the page
Chrome, Firefox, Safari, desktop and one phone, in a worker. Load and parse time, compute time, peak memory, native/wasm ratio, bytes fetched; real OOM n (wasm32 ceiling). Parity with native at every n. **Design decision made now:** the page fetches precomputed condensed matrices and ordinations per species (about 50 MB at 5k in f32); sketches (about 3 MB per sample in a bundle, about 15 GB at 5,000) are loaded only for a query or added samples. Pass: time to first plot.

### T6. Cross-species behaviour
Pool the 7 species plus distant controls: do 1 - ANI and containment separate within-species, between-species and between-genus, and where do they saturate. Pass: a written rule for when two samples are comparable (expected: group by species before ordinating).

### T7. Query tier
Hash-to-sample inverted index for 2,002 and replicated 5,000 samples: size, build time, "closest samples to this sketch" and "contains this sequence" latency versus brute force; sharded by hash prefix (no Range on the host). Pass: under 2 s at 5,000 with a per-species index small enough to hold. Waits on T6.

### T7b. Browser-side database
Whole-file `fetch` from S3 (CORS only, no Range) into DuckDB-wasm (`registerFileBuffer`, Parquet catalog) and into Rust/WASM (sorted delta-coded inverted index shard), thrown away on tab close. 300 and 1,000 samples per species: memory, load time, query latency, desktop and phone. About 90k hashes per *C. neoformans* sample means 27M rows at 300 samples; cross-species tables are too large, so shard per species or hash prefix. Expectation, not measured: DuckDB for metadata, Rust for the hash index.

### T8. Hosting and transfer
NRP S3 `s3://stajichlab/sketch/fungi/` (owner's `~/.s3cfg`, laptop only, `s3cmd -P`), throwaway file first: CORS from `https://gmg.fungalgenomes.org`, exposed `ETag`/`Content-Length`, Range support, concurrent GET throughput, egress limits. Pass: first view of a precomputed species under 30 s on a good connection (met by T5's design, not by fetching sketches).

### T9. Hypothesis tests
*(review)* No PERMANOVA exists in the code yet. Build it with PERMDISP (unequal spread is mistaken for a location effect when groups are unbalanced, as in *P. nodorum*). Permute within study (studies are clustered and often clonal, so samples are not exchangeable), not all labels; add planted-effect sets to measure power. Pass: calibrated p-values and a confounding warning.

### T10. Heterozygosity, ploidy and contamination at the capped depth *(review, new)*
At 20-30x, heterozygous k-mers sit near 10-15x, close to the abundance cutoff, and the shoulder test needs about 12x k-mer coverage. Mix reads in silico: a contaminant at 1/3/10/30% (a second species), and a diploid such as *C. albicans*, at each depth. Report both distances (the tree uses mash, PCoA uses containment ANI, so they disagree under extra content) and flag recall (`size_high`, `abundance_shoulder`, `unrelated_to_cohort`). May change go/no-go for the depth cap.

### T11. Reads versus assemblies *(review, new)*
The database will mix SRA read sketches and assemblies. Take strains with both (BFD/Cocci sets) and measure the offset in both distances; decide whether they may share an ordination.

### T1. Acquisition cost (cut down, last)
20 runs x 2 methods (current `parallel-fastq-dump` + sampler; `fasterq-dump`), stratified by run size, timed at several hours of the day. Add ENA streaming truncation only if T2 shows first-N reads are acceptable. Record wall, CPU, bytes in, scratch peak (fasterq-dump scratch is about 10x the .sra size), failure rate. A 300 Mb genome at 30x needs about 9 Gb kept, several times the cost of the small ones.

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
