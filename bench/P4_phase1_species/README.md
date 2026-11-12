# P4: phase 1 of the PopSeq panel plan, one full species (*Brettanomyces bruxellensis*, taxid 5007)

All 300 panel isolates of one species (299 run, 1 deferred), fetched at 25x by stride and sketched, a few fetches at a
time per node, on HPCC, 2026-11-11 (13:00 to 14:57 UTC). Plan: `docs/plan_popsra_scaleup.md`, phase 1.
`summarize.py` rebuilds `results/summary.md` from the raw files in `results/`; every number below is in it.

## The gate (from the plan)

| Gate | Result |
|---|---|
| At least 95% of runs `ok` | **299 of 299 `ok`, exit code 0 for every fetch**; every run verified by the full check (stride), no fallback |
| Share routed to ENA | 2 of 299 (0.7%); 297 by AWS. In the 10-isolate *C. auris* pilot it was 2 of 10 |
| Share failing the sampling check | 0 of 299 (but see finding 5) |
| Transfer rate per fetch | AWS median 18.5 to 20.0 MB/s (min 1.7, max 30.2); ENA 0.62 and 0.81 MB/s (n = 2) |
| Per-run cost | median 449 s per fetch (min 47, p90 625, max 4,362); 40.9 task-hours in all; peak memory per job 5.8 to 20.0 GiB |

198 GB were transferred for 296 Gb of runs (`.sra` median 0.42 bytes per base over 297 runs), of which 98 Gb of bases
were kept. 299 sketches of 3.1 to 9.5 MB (median 5.1 MB); no reads were kept.

## Findings

1. **The pipeline works at this scale.** 299 of 299, with the node packing described below.
2. **My cost model was low, and the time is not where I thought.** The median fetch took 0.74 s per Mb of run
   (estimate from three deep runs: 0.19 to 0.40), and 40.9 task-hours in all against my 16 to 33. The transfer is
   a small part: 21 s of a 449 s median fetch, conversion 50 s, the sketch 13 s. **80% of each fetch (median 326 s) is time
   the sidecar timers do not name**: the plan lookups, the read counts of the full files, the evenly spaced cut and the checks
   after it. They were not timed separately, so which dominates is UNMEASURED, but it is the place to look to make stride cheaper.
3. **Memory is tight at the top.** Task 0 peaked at 20.0 GiB of its 20 GiB request (99.98%); it held the two largest runs
   (7.7 and 7.6 Gb). The other nine jobs peaked at 5.8 to 10.7 GiB. It finished, but a 20 GiB request is not safe for chunks with large
   runs: the next phase needs a larger request, or a size cap per chunk. Which component uses the memory is not known.
4. **Packing did not slow the transfer, probably because transfers rarely overlap.** The median AWS rate was 18.5 MB/s with
   8 fetches on a node (r12), 20.0 with 8 (r39) and 18.9 with 16 (r11). But a transfer is only about 5% of a fetch's time, so
   few fetches were transferring at the same moment on any node. This does not contradict PopSeq's finding (speed per
   fetch depends on how many fetches share a node); it is a different regime. The rates have a long low tail (min 1.7 and 1.9 MB/s).
5. **The sampling check passed every sample, but with no margin on one.** Size ratio 0.50 to 2.35 (median 1.09), coverage ratio
   0.04 to 0.92. The smallest size ratio, 0.50, passes a threshold of 0.5 by a hair. Under stride the check is not needed
   (it exists to catch the bias of contiguous windows), but if windows with the fallback were used, genuine variation in
   this species (finding 6) would sit right at the threshold: the thresholds, from one species, do not transfer to this one.
6. **42% of the samples are flagged by the audit, and the data suggest the catalog's genome size is too small for many of them.**
   125 of 299 are flagged: `low_coverage` 107, `size_high` 16, `unrelated_to_cohort` 6, `size_low` 2. They are concentrated:
   84 of the 169 isolates of PRJEB71961 are flagged. The `low_coverage` samples have a median k-mer coverage of 2.1 against 14.0 for
   the unflagged, although every one received 25.0x of the catalog's 13.15 Mb. The k-mer histograms
   (`results/kmer_hist_sample.txt`, six flagged and four unflagged samples) show why: the flagged samples have a broad, shallow hump
   at abundance about 4 to 9, 74 to 125 thousand singleton k-mers and 172 to 213 thousand distinct k-mers, against a peak at 16 or
   more, 23 to 28 thousand singletons and about 100 thousand distinct k-mers in three of the four unflagged samples shown (the fourth,
   SRR10551669, has 120 thousand singletons and 197 thousand distinct k-mers, and a peak above 16). That is roughly a third
   of the per-k-mer coverage and twice the distinct k-mers from the same 25x of 13.15 Mb, in ten samples. **That fits genomes larger or more
   heterozygous than the reference (polyploid or hybrid strains), and I cannot separate it from a higher sequencing error rate with
   these data; it is a hypothesis.** Also, the audit's coverage estimate of 2.1 is lower than the hump visible in the histograms
   (about 4 to 9) for the six flagged samples shown; why is not determined.
7. **Collapse and label tests.** `--collapse mean --collapse-within study`: 299 samples to 122 representatives (largest cluster 68),
   noise floor 0.0055. Only `country` could be tested (n = 66 with a country, 10 groups): R2 0.055, p = 0.835, stratified by study.
   One species and mostly one study make this an example, not a result about the species.

## Limits

- One species, and one with variable ploidy: the cost, memory and flag numbers will differ for other species.
- The node packing was uneven (8, 16 and 8 fetches on three nodes), and not designed as an experiment.
- The 125 flagged samples are kept in the sketch set and marked; I have not decided whether to exclude them from an analysis.
- `SRR3927079` (47.9 Gb) is deferred: stride on it would take hours (`manifest_deferred.csv`).

## What it means for the next phase (not decided)

Raise the memory request or cap the run size per chunk; update the cost model to the measured 0.74 s per Mb (median) and look at the
unnamed 80%; decide how flagged samples are handled and whether the genome size used for the coverage target should come from the
sample's own k-mer content; do not rely on the sampling check's thresholds for polyploid species.

## Files

`manifest.csv`, `manifest_deferred.csv`, `summarize.py`, `results/` (`summary.md`, `results.tsv`, `exits.tsv`, `sampling_check.tsv`,
`hook_times.tsv`, `sacct.txt`, `downloads_p4.csv`, `sidecars.tgz`, `labels.tsv`, `kmer_hist.py`, `kmer_hist_sample.txt`, and the `popsketch run`
outputs `popsketch.json`, `qc.tsv`, `collapse.tsv`, `collapse_labels.tsv`, `group_tests.tsv`, `split_half.tsv`, `clusters.tsv`).
On HPCC the work is in `~/projects/GmG/bench/P4_phase1_build/bench/P4_phase1_species/work` (purgeable once this is merged).
