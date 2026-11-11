# Plan: sketching the PopSeq population panel (412 species, 34,673 isolates)

Status: draft, 2026-11-11. It says what is measured, what is arithmetic from measurements, and
what is UNMEASURED. It outlines how to scale up when more compute is available; it does not
start anything. Decisions for the owner are listed at the end.

## The input

PopSeq panel `panels/pop_sra_ge10_cap300/` (PopSeq PR #39, issue #36): species with at least 10
paired Illumina isolates, at most 300 isolates each, one run per isolate, ENA URLs in
`samples.csv`, depth and labels in `metadata.csv`. Each species is its own population: no pooled
ordination. Checked from the panel files (2026-11-11):

| Quantity | Value |
|---|---|
| Species / isolates | 412 / 34,673 |
| Isolates per species | median 36; 47 species have the full 300 |
| Whole runs, all isolates | 135.3 Tbases |
| Kept at min(25x, whole run), isolates with a genome size | 33.8 Tbases (32,785 isolates) |
| Isolates with no genome size (36 species) | 1,888: the coverage tier cannot target them |
| Runs shorter than 25x (whole run kept) | 2,843 of 32,785 |
| Isolates below 10x (kept in the panel) | 549 |

## What is measured

| Fact | Value | Source |
|---|---|---|
| `popseq fetch` at every tier, one node | 11 of 11 acceptance steps passed; coverage target reached to 3 decimals | PopSeq PR #38, `benchmarks/2026-11-09_fetch_acceptance_*.tsv` |
| AWS `.sra` rate per task, HPCC | 3.0 to 7.6 MB/s in 8 acceptance rows (one node, one afternoon, unknown load); 10.3 MB/s median over 4,211 MATPredict runs | PopSeq spec 002 |
| ENA rate per task, HPCC | 0.97 MB/s | acceptance row, ERR1069479 |
| Laptop pilot (below) | ENA 0.22 to 0.30 MB/s per task at 3 tasks; AWS 2.7 MB/s by `curl`, but the laptop has no SRA toolkit | `bench/P1_pilot_coverage/` |
| gzip FASTQ bytes per base | median 0.55, range 0.34 to 0.83 (n = 17 runs). The spec's divisor (bases / 3.5) predicts 0.29: **the divisor is low by 1.9 times, 1.2 to 2.9 times** | benchmark runs and pilot plan, PopSeq `popseq/benchmarks`, `bench/P1_pilot_coverage/plan.csv` |
| `.sra` bytes per base | 0.28, 0.49, 0.65 (n = 3 only) | acceptance rows |
| Sketch size | 4.5 MB whole, 11.4 MB with split halves, per sample | `bench/A1_measure`, #167 |
| `popsketch run` on one species of 300 | about 2 GB resident, 7 s | A1 |
| Distance method at this size | the fast bitset path, with a large margin (0.59 GiB of 4 GiB for 1,822 pooled samples) | D1b, #169 |
| Noise model | sampling noise is 10 to 100 times below the split-half floors; scaled 1000 gives sketches 5.0 times smaller | A2, #177 |

## The laptop pilot (this session)

Three species, one run each, `popseq fetch --tier coverage --coverage 25`, from `$SCRATCH`, ENA
only (the laptop has no SRA toolkit), three tasks at once (declared concurrency 3).
Script and results: `bench/P1_pilot_coverage/` (`select_pilot.py`, `run_pilot.sh`, `results.tsv`).

| Run | Species | Result | Coverage | Transfer |
|---|---|---|---|---|
| ERR12524156 | *Brettanomyces bruxellensis* | ok, verified | 25.000x of 13.15 Mb | 293 MB, 1,022 s, 0.287 MB/s |
| SRR31272388 | *Candida albicans* | ok, verified | 24.999x of 14.28 Mb | 329 MB, 1,099 s, 0.299 MB/s |
| ERR4669779 | *Nakaseomyces glabratus* | `verify_failed` (exit 12) on the first try; the retry hung and was stopped (exit 13) | none | first try: mate 2 arrived as 797 bytes, not 154,936,766. Retry: mate 1 stopped growing at 83.6 of 152 MB for over 3 minutes |

What the pilot shows, and what it does not:

- **The pipeline works end to end on a laptop**, from the panel's ENA URLs to a verified,
  coverage-targeted pair of FASTQ files and a sidecar.
- **The coverage tier does not reduce the transfer.** The whole run is downloaded (293 MB) and
  202 MB are kept (69%, as 25x of a 35x run). It reduces what is stored and what is sketched.
- **Verification caught a bad transfer.** One mate arrived as a 797-byte stub; the MD5 check
  failed it and nothing was kept. The outcome is `verify_failed`, but a size that is not ENA's
  `fastq_bytes` is a short transfer, which is transient (exit 13), and a size check before the
  MD5 would say so. A follow-up for PopSeq; see below.
- **A transfer can hang with no timeout.** The retry's `curl` was alive for 29 minutes and its file
  stopped growing; a second, unrelated probe to ENA also hung. `popseq fetch` calls `curl` with
  `--retry 3` but no stall limit, so a stalled connection waits until the server closes it, and a
  restart begins again from byte 0 (the ENA transfer does not use `curl -C -`, although the spec
  says it resumes where it can). I stopped it by hand; it ended cleanly (exit 13, scratch removed,
  nothing partial kept). The cause of the stall is not known: this laptop's route to ENA, or ENA.
- **A laptop is not a platform for this.** At the measured 0.3 MB/s per task and 3 tasks,
  74 TB of FASTQ is about 2.6 years. The pilot is a test of the pipeline, not of throughput.
- Two of three runs finished, so this is n = 2 for the successful rates (0.287 and 0.299 MB/s) and n = 1 for each failure.
- Not shown: sketching (the pilot did not run the hook), the SRA routes (no toolkit here), any
  rate under a realistic cluster load.

## Arithmetic for the full panel (assumptions stated; not measurements)

Transfer volume is the **whole runs**, not the 33.8 Tbases kept: `fastq-dump -N/-X` needs the
whole `.sra`, and ENA serves whole files.

| Route | Volume | Basis |
|---|---|---|
| ENA FASTQ | 74 TB (46 to 112) | 135.3 Tbases x 0.55 (0.34 to 0.83) bytes per base, n = 17 |
| AWS `.sra` | about 66 TB (38 to 88) | 135.3 Tbases x 0.49 (0.28 to 0.65), **n = 3 only** |

Wall time for the transfer alone is `volume / (tasks x rate per task)`:

| Scenario | Wall time |
|---|---|
| AWS, 20 tasks at 5 MB/s (the documented starting cap, UNMEASURED) | 7.6 days |
| AWS, 100 tasks at 5 MB/s | 1.5 days |
| AWS, 100 tasks at 10 MB/s | 0.8 days |
| ENA, 3 tasks at 0.3 MB/s (this laptop) | 957 days |

Whether 100 tasks hold 5 to 10 MB/s each is exactly what issue #27 asks; nobody has measured it.
Sketching time per run is UNMEASURED (the first number comes from phase 1 below). Storage if no
reads are kept: 34,673 x 11.4 MB = 395 GB of sketches (156 GB whole sketches only), at most 3.4 GB
per species.

## Plan, in phases with a gate each

**Phase 0: laptop pilot (done above).** Gate met: fetch works end to end; one transfer failure
was caught by verification.

**Phase 1: one species end to end on HPCC.** A species with 300 isolates and a genome size (for
example *C. auris*, 12.25 Mb). `popseq fetch --tier coverage --coverage 25 --exec <sketch hook>`
per run in a Slurm array, then `popsketch run --from-sketches` with `--collapse --collapse-within
study`. Needs: a small hook script that calls `popsketch sketch --id RUN` on the two FASTQ paths
(to write, in `scripts/`, with a test); split halves for every sample (`--floor-samples 1000`),
because the collapse noise comes from them. Measures: wall time, memory and sketch time per run
(T1), the failure rate by exit code, and the real bytes per base for this species. Gate: at least
95% of runs `ok`, and the per-run cost is known.

**Phase 2: about 10 species, 3,000 runs, ramping the concurrency.** The same array, with the
number of concurrent tasks raised 5, 10, 20, 40 while `benchmarks/downloads.csv` records MB/s and
outcomes with the declared `concurrency`. This is the measurement for #27, from real rows.
Gate: a concurrency cap chosen from rates and failure rates, not from a guess.

**Phase 3: the full panel in waves by species.** Smallest species first (median 36 isolates), then
the 47 species of 300. One manifest per wave; `popseq fetch` resumes from sidecars, so a failed
wave is rerun, not restarted. Retry by exit code: 13 with back-off; 12 once, then set aside;
10 and 11 recorded and skipped; 14 reported. The 36 species with no genome size wait for the
owner's decision (below).

**Phase 4: results and hosting.** One result folder per species (matrices, PCoA, PCA, tree,
clusters, label tests, `collapse.tsv`), a species index page, precomputed matrices for the
browser (the page cannot hold 300 sketches), and S3 hosting under `s3://stajichlab/sketch/fungi/`
(CORS for `https://gmg.fungalgenomes.org` is set; a real-browser check is still open).

## What more compute changes

- **Throughput is the constraint, not CPU.** One species fits a laptop's memory; the transfer of
  66 to 74 TB does not fit a laptop's network. More compute means more concurrent tasks on a
  cluster with fast egress, and the cap measured in phase 2.
- **Route.** AWS where NCBI has the run, ENA where it does not or holds it unpaired. ENA at about
  1 MB/s per task is the slow tail, so the fraction of runs that must come from ENA decides the
  tail of the wall time. UNMEASURED for this panel; phase 1 gives it per species.
- **Storage stays small** if reads are not kept (`--exec`). Keeping 25x of FASTQ would be about
  34 Tbases x 0.55 = 19 TB, which is why the hook matters.

## Gaps and risks

- **The spec's storage divisor is wrong** (0.55 observed against 0.29 predicted). Any storage
  figure that used it is low by about two.
- **Short transfers are reported as `verify_failed`** (found in the pilot). A size check against
  ENA `fastq_bytes` before the MD5 would classify them as transient.
- **A hung transfer has no timeout and does not resume** (found in the pilot): add a stall limit
  (`curl --speed-limit/--speed-time`) and `-C -`, so a stall ends as transient and a retry continues.
- **`popseq fetch` has no preflight for missing tools.** On a host without the SRA toolkit it
  would download a whole `.sra` and then fail at conversion; the pilot avoided this by giving a
  plan without NCBI counts. A tool check would skip those routes without a transfer.
- **Collapse needs split halves** for every sample, or the noise falls back to the species median.
- **549 isolates are below 10x** and 2,843 runs are shorter than 25x: the panel keeps them; filter
  on `est_depth_x` if the sketch needs depth.

## Decisions for the owner

1. **The 36 species with no genome size (1,888 isolates).** Supply sizes, use the `windows` or
   `stride` tier for them, or leave them out of the first pass. I would not guess a size.
2. **The depth target.** 25x is the working choice; 20 to 30x was the earlier range.
3. **Where phase 1 and 2 run, and how many concurrent tasks may be used.** The cap is a policy
   decision as much as a measurement.
4. **Whether to file the four PopSeq follow-ups** listed under Gaps (size check before the MD5,
   stall timeout and resume, tool preflight, the storage divisor in spec 002 section "Size
   estimation"). The code is PopSeq's, so the fixes belong there.

## Updates

**2026-11-11: sampling and the per-node packing of phase 2.**

- **The coverage tier now samples by stride by default.** The phase 1 pilot (`bench/P2_phase1_pilot/`) showed that contiguous
  windows reach only `25 / depth` of the genome on a position-sorted run (two of ten runs came out as a 2.2 and a 2.3 Mb
  genome). PopSeq's `popseq fetch` (branch `fetch-stride`) keeps evenly spaced reads of the verified full run, with windows
  opt-in and a fallback from windows to stride when the sketch hook's check fails (`popsketch check`, exit 75). Tested on
  the biased runs in `bench/P3_stride_fallback/`: the genome estimate went from 2.2 and 2.3 Mb to 11.9 and 12.2 Mb, and the
  fallback did not download again. n = 3 runs of one species.
- **It costs more.** Stride took 3.5 to 4.0 times the seconds of windows on the same three runs, and a peak of 3.1 GB. Phase
  2's compute estimates must use the stride cost, not the phase 1 windows cost. A sorted run taken by windows with the
  fallback costs about what stride costs; a run in random order costs the windows time. Which is cheaper overall depends on
  the share of sorted runs: 2 of 10 in the one species measured, unknown for the panel. Measure it in phase 2.
- **The thresholds of the check are untested beyond one species** (size under 0.5 of the genome, k-mer coverage over 2.5
  times the expected). Phase 2 should log the check for every run, so the thresholds can be examined across species.
- **Pack the fetches per node** (PopSeq issue #27, spec 002 "Concurrency sweep results"): speed per fetch depends on how many
  fetches share a node, about 17 MB/s at 2 to 5 per node and about 5 MB/s at 20 or more; Slurm packed 52 of 80 array tasks
  onto one node, and the `%N` array limit does not control that. For phase 2, limit the tasks per node (for example a fixed
  number of parallel fetches inside each job). That arrangement is unmeasured.
