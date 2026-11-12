# P6: where do time and memory go in a chunk job? (2026-11-12)

Re-ran two chunks of 2 *B. bruxellensis* runs (`manifest.csv`) with `popseq fetch --profile` (PopSeq #45) through the chunk script
(GmG #186): chunk 0 = the two largest runs of phase 1's task 0 (the one that peaked at 20.0 of 20 GiB), chunk 1 = two smaller runs.
PAR=2, 8 threads each, 16 CPUs, 32 GB requested (so the peak could show), stride sampling, coverage 25. Job 32024792. Reproduce with `run.sh`.
`results/profile.tsv` is the per-fetch profile, `results/sacct.txt` the Slurm accounting.

| Run | Gb | .sra downloaded | s | transfer | convert | cut | hook | tree peak MiB (stage) |
|---|---|---|---|---|---|---|---|---|
| ERR2004534 | 1.5 | 1.05 GB | 355 | 48 | 72 | 200 | 11 | 675 (convert) |
| SRR945538 | 1.1 | 0.74 GB | 369 | 37 | 56 | 238 | 11 | 646 (convert) |
| ERR1102659 | 7.7 | 3.46 GB | 1,389 | 171 | 415 | 676 | 16 | 595 (convert) |
| SRR065689 | 7.6 | **51.46 GB** | 3,984 | **2,833** | 367 | 641 | 13 | 535 (convert) |

(.sra sizes are from phase 1's `downloads.csv`, same runs and route. Every run passed its check.)

## Findings

1. **The stride cut is the biggest named stage** (200 to 676 s, about half of a fetch that does not wait on a transfer). This is the
   "80% unnamed" of phase 1. `other_s` is 0.5 to 2.9 s.
2. **Process memory is small.** Whole-tree peak 535 to 675 MiB, all in `convert` (parallel-fastq-dump); Python 47 to 126 MiB.
   `parallel-fastq-dump` is the largest stage but it is not the memory problem. (Sampled every 2 s: a lower bound.)
3. **The 20 GiB (phase 1) and 33.5 GB (here, at the 32 GB cap) are not process memory.** sacct MaxRSS of the batch step was 33,550,564 K for
   chunk 0 and 5,488,664 K for chunk 1, against tree peaks of 0.5 to 0.7 GiB. It follows the size of what the job reads and writes
   (page cache), as the PopSeq session found for its own runs. So memory requests need not be sized from MaxRSS.
4. **SRR065689's object at AWS is 51.5 GB for 7.6 Gb of bases** (15x the usual 0.4 bytes per base): 2,833 s of transfer, 71% of its fetch,
   and 51 GB written to `$SCRATCH`. Why (probably a different `.sra` flavour with original-format data) is not established.
   It is the real outlier of this chunk, for time and for `$SCRATCH` space.

## Consequences (for the owner to weigh)
- Request memory for the process (about 1 GiB per fetch measured here, plus the sketch hook, which is inside these numbers), not from MaxRSS.
  This is for 4 runs of one species: UNMEASURED beyond that.
- A per-run `.sra` size check before transfer (the PopSeq follow-up already proposed) would catch a 51 GB object before it costs 47 min.
- Only these runs were profiled: whether the cut share (about 50%) holds across species is UNMEASURED.
