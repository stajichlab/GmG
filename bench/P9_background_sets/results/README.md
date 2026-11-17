# P9 results: background sets (2026-11-14 to 2026-11-16)

The *Coccidioides* pair (2 x 50 runs) and the *A. fumigatus* pool (150 runs) from `analysis/popsets`, fetched with `popseq fetch --profile` (PopSeq 0.2.0, coverage 25, `--sampling auto`, partition `stajichlab`, node r11), 8 chunks of up to 32 runs plus one retry pass. Files: `profile.tsv` (every attempt, per fetch), `results.tsv`, `exits.tsv` (every attempt), `sampling_check.tsv`, `hook_times.tsv`. Summary made with `../summarize.py`.

## Outcome
**250 of 250 runs ok** after the retry pass (50 + 50 + 150 per set). The first pass had transient ENA failures (exit 13); `exits.tsv` keeps every attempt, `summarize.py` reads the latest one per run. `manifest_deferred.csv` holds only its header: no run is deferred in the file (an earlier note of "one 48.5 Gb run deferred" does not match it).

## Cost of the 250 ok fetches (median, max)
| Stage | seconds |
|---|---|
| whole fetch | 365, 5,096 |
| transfer | 89, 4,132 |
| convert (parallel-fastq-dump) | 226, 291 |
| count | 12, 63 |
| cut | 0.2, 855 |
| hook (sketch + check) | 28, 31 |

Sampling: **windows for 249 runs, stride for 1** (the cut costs 855 s for that one run, the stride path). Process tree peak memory: median 149 MiB, max 455 MiB. The slowest fetches are ENA transfers (4,132 s max), not conversion or sketching.

## Sets sketched (fixed popsketch 0.3.2, `run --from-sketches`, 1 per run)
| Set | in analysis | median size (Mb) | gap verdict | flags |
|---|---|---|---|---|
| *C. immitis* (50) | 49 / 50 | 26.2 | **clear cut into 2** (jump 10.98 floors) | 1 unrelated_to_cohort |
| *C. posadasii* (50) | 49 / 50 | 26.0 | no clear cut (jump 0.22) | 1 size_high |
| *A. fumigatus* pool (150) | 149 / 150 | 28.0 | no clear cut (jump 1.18) | 1 low_coverage |

The *Coccidioides* results match the earlier 559-sample analysis (*C. immitis* one cut into 2; *C. posadasii* none). Not shown: whether the two *C. immitis* clusters follow geography, and the *A. fumigatus* pool against labels or clades. Sketches and `popsketch.json` are on the HPCC under `bench/P9_background_sets/out/{work/sketches,sets_out}/` (purgeable).
