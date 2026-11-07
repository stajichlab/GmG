# Plan: a test framework for "popsketch at scale" (draft v1, 2026-11-07)

## The question

Can we, for thousands of strains across several species: download each run, subsample to about 20-30x, sketch it at a good size/sensitivity trade-off, and then give a browser (WASM) page the means to make tables, PCoA/PCA, trees, hypothesis tests and a query over the sketches? This plan is not the build plan (see `plan_web_sketches.md`). It is the set of measurements that tells us (1) what is possible and how well, and (2) what resources it takes, before we commit.

## Corrections to what we believed

- The fetch script (`03_fetch_sketch.sbatch`) keeps every Nth read pair of the whole run (systematic subsample), not the first reads. It must download and decompress the whole run, so cost scales with run size, not with the 1 Gb we keep.
- "20-30x" is a per-species base cap (genome size x 25), not a constant. The fixed 1 Gb used so far gives 20-100x for 18-35 Mb genomes and would give about 3x for a 300 Mb genome.
- Known and measured: 8-13 min wall on 4 cores per run at 1 Gb (including download); native analysis 1 min at 1,500 samples and 4.5 min at 3,000 (synthetic hashes); 8 MB per sample now (about 1 MB planned compact). Everything else below is an estimate until a test replaces it.

## Method

Each test has: a question, a dataset, a metric, a pass criterion set **before** running, and a resource record. Tests are ordered cheap to expensive, and a failed gate stops the later tests that depend on it.

**Resource record.** Every job writes one line to `analysis/bench/results.tsv`: test id, sample id, wall s, CPU s, max RSS, bytes downloaded, bytes written, scratch peak, node type, exit state. Slurm `sacct` fields plus `/usr/bin/time -v` for anything run outside Slurm; browser runs record `performance.now`, `performance.memory` where available, and transferred bytes. Costs are extrapolated from fits (cost = a + b x run size), not from single points.

**Data.** Existing: the 7 SRA population sets (2,002 sketches, 1 Gb each), the Coccidioides and B. dendrobatidis sets with known problems, the 12-run deep sets of C. neoformans and Z. tritici. Needed: runs of a larger genome (100-300 Mb, species to be chosen) so the depth cap is tested where it bites. Synthetic data only for pure scaling (n), never for accuracy claims.

**Reference.** Never a very deep whole run (shown unreliable). Use the 6 Gb subsample for small genomes and report the dependence of Pearson r on the span of diversity in the set, not r alone.

## Tests

### T1. Acquisition cost and failure rate (resources)
60 runs stratified by size (0.5, 2, 8, 30 Gb). Compare: current `parallel-fastq-dump` + awk subsample; `fasterq-dump`; ENA FASTQ over HTTPS streamed through `head`-style truncation (cheap, but biased: SRA FASTQ is position-sorted for old runs); `prefetch` then dump. Record wall, CPU, bytes in, scratch peak, failure rate, retries. Output: cost model per method, and whether streaming a truncated prefix is ever acceptable (T2 decides).

### T2. Subsample design and depth (what is possible)
From one download per run produce: every-Nth at 5, 10, 15, 20, 30, 50x; first-N reads; a random subset. 12 runs each of C. neoformans, Z. tritici and one genome >= 100 Mb. Metrics against the reference: Pearson r and mean absolute difference of the pair distances, nearest-neighbour agreement, ARI of the clusters, estimator "peak found" rate, noise floor, sketch size. Pass: nearest-neighbour agreement >= 95% and mean difference below the noise floor at the chosen depth; the smallest depth that passes per species gives the per-species cap rule (cap = f x genome size) and the saving in T1 cost.

### T3. Sketch parameters: size versus sensitivity
Same reads, k in {17, 21, 25}, scaled in {100, 200, 500, 1000}. Metrics: bytes per sample (hash count x 5 B), distance fidelity versus k=21/scaled 100, and the resolution: distance divided by split-half noise floor for known close pairs, and recovery of known clades (Cocci, Bd sets). Pass: the cheapest setting that keeps the cluster and NN results of T2 within tolerance. Also gives the real size for the hosting plan.

### T4. Native scaling in n
Real sketches replicated by resampling to n = 500, 1k, 2k, 5k, 10k, plus real 2,002. Time and RSS for distances, tree, hash PCA, PCoA, clustering. Check the n^2 matrix memory (10k samples = 800 MB as f64) and where the eigen-solver dominates; evaluate randomized/partial eigen-solvers and a condensed (f32, upper triangle) matrix if superlinear. Pass: a stated n_max for a 4 GB / 60 s budget.

### T5. WASM and the page
Same n in Chrome, Firefox and Safari, desktop and one phone; in a worker. Record load/parse time, compute time, peak memory, native-to-wasm ratio, and bytes fetched. wasm32 has a 4 GB address-space ceiling, so find the real out-of-memory n. Parity with native on every n (extends the existing wasm parity test). Pass: n_max_browser per device class, and time to first plot.

### T6. Cross-species behaviour
Pool all 7 species (+ distant controls): do 1 - ANI and containment separate within-species, between-species and between-genus, and where does the ANI estimate saturate at scaled 200? Pass: a written rule for when two samples are comparable, and whether the page must group by species before ordinating (expected: yes).

### T7. Query tier
Build a hash-to-sample inverted index for the 2,002 sketches and for a replicated 5,000. Index size and build time, latency for "closest samples to this sketch" and "samples containing this sequence", versus brute force over sketches. Shards by hash prefix so it works on a host without byte-range requests. Pass: latency under 2 s at 5,000 samples with an index small enough to hold per species.

### T8. Hosting and transfer
On NRP S3 (`s3://stajichlab/sketch/fungi/`, owner's credentials, laptop only): CORS from `https://gmg.fungalgenomes.org`, `ETag`/`Content-Length` exposed, whether Range works, concurrent GET throughput from a browser, time to first view for 100/1,000 samples, egress limits. Use a throwaway file before any real upload. Pass: first view of a species under 30 s on a good connection.

### T9. Hypothesis tests
PERMANOVA (and nearest-neighbour label purity) on the distance matrix by label (country, study, host) on sets with known structure and on label-permuted copies. Runtime in wasm, calibration (false positive rate near the nominal level under permutation), and sensitivity to the sampling imbalance already seen (P. nodorum: 187 of 265 from Australia). Pass: calibrated p-values; warnings when groups are confounded.

## Resource ledger (to fill in)

| item | per sample | per 5,000 samples | source |
|---|---|---|---|
| download+subsample wall, core-hours | T1 | | |
| scratch peak | T1 | | |
| sketch bytes, compact | T3 | | |
| analysis time native / wasm | T4 / T5 | | |
| browser memory | T5 | | |
| storage + egress | T8 | | |
| index | T7 | | |

Planning estimates (not measured): T1+T2 about 250 core-hours, T3 reuses T2 reads, T4 about 10 core-hours, T5-T9 laptop. Under 1 TB transient scratch, about one week wall, run when HPCC load is low.

## Decisions the results feed

- depth cap rule per species, and the subsample method;
- k, scaled and file format (size per sample);
- per-species n_max for native and browser, and whether a server tier is needed beyond that;
- whether the query tier is a static sharded index or a service;
- host, and cost per 1,000 views.

## Order and gates

1. Harness and resource record (`analysis/bench/`), T8 on a throwaway file, T4 (no HPCC).
2. T1, then T2 and T3 from the same downloads (HPCC, low load).
3. T5, T6, T9 on the results; T7.
4. Write the ledger, decide the build plan.
A gate fails if a pass criterion above is missed; that changes the plan, it is not a reason to relax the criterion afterwards.

## Risks

- Reference bias (see Method); small samples of runs per species; one large genome may not represent others.
- NCBI/ENA throughput varies by hour; time T1 at several hours of the day.
- Browser memory differs a lot by device; the phone result may force a smaller n.
- Cross-species result may show the sketch cannot compare beyond genus, which limits the "multiple species" view to a catalog and search.
