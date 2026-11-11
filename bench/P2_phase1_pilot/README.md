# P2: phase 1 mini-pilot of the PopSeq panel plan (10 isolates of *Candidozyma auris*, HPCC)

`popseq fetch --tier coverage --coverage 25 --require-pairs --exec scripts/popsra_sketch_hook.sh`, one run per Slurm
array task (`scripts/popsra_fetch_array.slurm`, 5 at once, declared concurrency 5), then `popsketch run
--from-sketches --collapse mean` over the 10 sketches. Plan: `docs/plan_popsra_scaleup.md`, phase 1. The runs are
the first run of each of 10 studies in the panel (`select_phase1.py`, depth 25 to 200x). Raw results are in
`results/`; `summarize.py` rebuilds the table below from them.

## Measured (one species, 10 runs, HPCC, 2026-11-11)

See `results/summary.md` for the full table. Summary:

| Quantity | Value |
|---|---|
| Outcome | 10 of 10 runs `ok` after one rerun (see below); reads never kept (0 FASTQ files left) |
| Coverage reached | 24.996x to 25.033x of 12.25 Mb |
| Transfer, AWS route (8 runs) | 4.0 to 9.7 MB/s per task, median 5.8, at 5 tasks at once |
| Transfer, ENA route (2 runs) | 0.86 and 0.92 MB/s per task |
| Runs that needed ENA | 2 of 10: one absent at NCBI (ERR17725644), one held unpaired at NCBI (ERR2299868) |
| Task wall time | 1.5 to 27.8 min, median 2.8; the two ENA tasks took 10.4 and 27.8 min |
| Peak memory (fetch + sketch) | 0.62 to 2.18 GB, median 1.41 (8 GB was requested; 4 GB would have done) |
| Sketch step | 10 to 17 s per run |
| Sketch files per run (whole sketch, both halves, metadata) | 1.3 to 14.3 MB, median 4.4 |

The transfer was 0.25 to 1.5 GB per run: 7.7 GB for the ten successful transfers, 8.6 GB with the rejected first
attempt of ERR2299868 (0.89 GB), for 3.1 Gb of bases kept. The whole run is downloaded whatever the tier keeps. n = 10 runs of one species: not a basis for the whole panel.

## Two things the pilot found

**1. `--require-pairs` is needed for this panel, and it was missing from my first array script.** ERR2299868 is
single-end at NCBI (`spots_with_mates = 0`), but `fastq-dump` returned mate files. `popseq fetch` refused them
(`verify_failed`, exit 12), correctly. With `--require-pairs` the run is taken from ENA as verified pairs
(rerun, job 31992298: 25.000x, 1.08 GB at 0.92 MB/s, 27.8 min). The array script now always passes it.

**2. Contiguous windows are biased on position-sorted runs.** Two of the 10 runs (ERR2300774 and ERR12543635) were
flagged `size_low` by the popsketch audit: an estimated genome size of 2.2 and 2.3 Mb against a cohort median of
12.1 Mb, and a k-mer coverage of 106 to 107x against 12 to 20x for the other eight. If a run is sorted by
genome position, 25x taken as contiguous windows covers only 25 / depth of the genome, here
25/134 x 12.25 = 2.3 Mb and 25/140 x 12.25 = 2.2 Mb, within 5% of the sketch's estimates. The other eight
runs estimate 11.9 to 12.4 Mb. (That the two runs are position-sorted is inferred from this match; the read order of
the `.sra` was not inspected.) The sketches of those two are incomplete (1.3 and 1.5 MB), and their distances to the
others are not reliable. Taking every Nth read (`stride`) does not have this problem, because it samples every region of
a sorted file; contiguous windows cannot, however many there are, because the genome fraction they reach is fixed by
the bases kept and the depth. This is a flaw in the coverage tier as I specified it (PopSeq spec 002, from the
MATPredict `windows` pattern), and the reason GmG's original sampler used every Nth pair. The audit caught both
cases; without it they would have entered the sketch set unnoticed.

## What the collapse showed

`--collapse mean`: 10 samples to 4 representatives, noise floor 0.00412 (the 90th percentile of 10 split-half
distances, with halves of 12x). *C. auris* is near-clonal within clades, so some merging is expected, but the two
biased sketches and the pessimistic floor (see plan D1c) make this count unreliable; it is not a result about the species.
Ten isolates from ten studies is also too few to test `--collapse-within`.

## Files

`select_phase1.py` and `manifest.csv` (the 10 runs); `summarize.py`; `results/` (`results.tsv` the fetch result lines,
`hook_times.tsv`, `sacct.tsv`, `downloads_phase1.csv` the log rows, and the `popsketch run` outputs `popsketch.json`,
`qc.tsv`, `collapse.tsv`, `clusters.tsv`, `split_half.tsv`, `ani.tsv`). On HPCC the work is in
`~/projects/GmG/bench/P2_phase1_build/` (purgeable once this is merged).
