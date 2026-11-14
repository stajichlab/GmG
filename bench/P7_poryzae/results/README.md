# P7 chunk 0: the first 32 *P. oryzae* runs (2026-11-12)

Job 32037075 (chunk 0 of `manifest.csv`: its first 32 rows, not a random sample), partition `stajichlab`, node r11, 16 CPUs, 12 GB
requested, 4 fetches at a time with 4 threads each, coverage 25, `--sampling auto`, PopSeq 0.2.0, `--profile`. Queued 30 min, ran 1 h 52 min.
Files: `chunk0/profile.tsv` (per fetch), `results.tsv`, `exits.tsv`, `sampling_check.tsv`, `hook_times.tsv`, `sacct.txt`.

## Outcomes
29 of 32 ok. All 29 were sampled by **windows** (`auto` found no sorted run in this chunk). Three did not finish:

| Run | Outcome | Why |
|---|---|---|
| ERR18066360 | transient (13), ENA | `curl: (18) transfer closed with 1.68 GB remaining` after 385 s |
| ERR2061624 | transient (13), ENA | `curl: (18) transfer closed with 2.84 GB remaining` after **3,720 s** (a stall of about an hour) |
| ERR3771233 | hook_failed (14) | sketch check refused it, also at stride: estimated size 16.0 Mb is 39% of the 41.0 Mb genome. The run holds only 354 Mb of bases (est. depth 8.6x) against a 1,025 Mb target, so the sample cannot cover the genome. Correct refusal, not a bug |

## Cost of the 29 ok fetches (median, max)
| Stage | seconds |
|---|---|
| whole fetch | 500, 791 |
| transfer | 83, 316 |
| convert (parallel-fastq-dump) | 337, 405 |
| count | 24, 39 |
| cut | **0.4**, 0.9 |
| hook (sketch + check) | 43, 72 |
| other | 2.2, 6.2 |

Memory: whole process tree peak median **200 MiB** (max 217; 23 runs in `convert`, 6 in `hook`); Python peak 64 to 136 MiB.
sacct MaxRSS of the batch step was 12.58 GB against the 12 GB request with no kill: page cache, as in P6. 12 GB was far more than needed.

## What it says
- With windows the cut costs 0.4 s, against 200 to 676 s with stride in P6 (B. bruxellensis). Conversion is now the largest stage (about two thirds).
- Despite about 8x the bases per run, a fetch (median 500 s) is close to phase 1's 449 s, and the chunk ran 1.9 task-hours for 32 runs
  (phase 1: 40.9 task-hours for 299; 10 chunks of this kind would be about 19). Extrapolation from 32 runs: UNMEASURED beyond chunk 0.
- The two ENA transfers failed in the middle and one stalled for an hour: the stall timeout plus `curl -C` resume follow-up for PopSeq is now evidenced.
