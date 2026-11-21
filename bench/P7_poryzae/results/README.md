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

# All 300 runs (chunks 0 to 9 plus one retry pass; 2026-11-16)

Files: `all/{profile,results,exits,sampling_check,hook_times}.tsv` (every attempt; `exits.tsv` keeps all, the latest per run counts). Summary made with `../summarize.py`.

## Outcome: 281 of 300 ok
| Latest outcome | Runs |
|---|---|
| ok (0) | 281 |
| transient (13) | 15 |
| hook_failed (14) | 4 (sketch check refused: run too shallow for the genome; correct refusals) |

The first pass left 17 transient and 3 runs with no exit line; the retry pass (20 runs) recovered 5. **All 15 that still fail died during the ENA transfer** after 2 s to 4,200 s (files of 2.6 to 6.9 GB), so more retries from the same host will likely fail the same way: this is the missing stall timeout and resume in PopSeq (stajichlab/PopSeq#51), not a GmG problem. Failed runs: ERR1448263, ERR1448265, ERR1448266, ERR1448268, ERR1448269, ERR18066360, ERR18066403, ERR18066418, ERR18066421, ERR18066425, ERR2061624, ERR2653642, ERR2653650, ERR2660592, ERR2662958.

## Cost of the 281 ok fetches (median, max)
| Stage | seconds |
|---|---|
| whole fetch | 565, 10,462 |
| transfer | 81, 8,850 |
| convert (parallel-fastq-dump) | 372, 862 |
| count | 24, 148 |
| cut | 0.4, 1,288 |
| hook (sketch + check) | 57, 90 |

Sampling: windows 265, **stride 14**, none recorded 2. Process tree peak memory: median 158 MiB, max 675 MiB (12 GB was far more than needed). The chunk 0 estimate of about 19 task-hours for the whole species is superseded by the measured run: see `sacct` of job 32116633 for the actual wall time (not tabulated here).

## The set, sketched (fixed popsketch 0.3.2, `run --from-sketches`)
281 sketches: 255 in analysis, median size 38.9 Mb (assembly 42 Mb), **no clear cut** (jump 0.86 floors, noise floor 0.0020), flags: abundance_shoulder 4, low_coverage 21, unrelated_to_cohort 9. Not shown: structure against geography or labels (the P7 manifest has study and country only); this is the candidate real population for stajichlab/GmG#155, not yet used for it.

# Update 2026-11-21: the 15 ENA failures re-fetched with PopSeq #55 and #56 (stall timeout, `curl -C -` resume)

All 15 runs listed above now finish, and the set is **296 of 300 ok** plus the 4 correct refusals (14). The re-fetch used PopSeq main `06539ed`,
`--sampling auto` (stride, coverage 25), 4 fetches at a time, species check `refuse` against `bfd_s10000.gmgi` (all 20 re-fetched runs: `match`; table
`all/species_match_retry.tsv`, one row per fetch). Files: the new rows are appended to `all/{results,exits,profile,hook_times,sampling_check}.tsv`.

- The retry manifest also held 5 runs that had already finished in the first pass (ERR2653636, ERR2653651, ERR2188724, ERR2660590, ERR2061621); they were fetched again, so
  those 5 have two or three rows in the files above; the latest row per run counts. This was a mistake of the retry list, not a finding.
- Pass 1 (4 h job limit) finished 6 of 20 (5 ok, ERR18066360 failed its md5) and hit the limit with the rest still downloading. Pass 2 (12 h limit, the 14 unfinished runs)
  finished all 14: 13 ok and ERR2188724 failed its md5. ENA served the 2 to 6 GB files at about 1 MB/s, cut connections (curl 18) and once sent an empty reply (curl 52);
  resume continued every time.
- **The two md5 failures (ERR18066360_2, ERR2188724_2) did not reproduce.** A plain `curl` of ERR18066360_2 matched ENA's md5; ERR2188724_2 matched after 12 and 7
  forced resumes in a stand-alone test; and the two runs passed on a third fetch (ERR2188724_2: 2 resumes, md5 ok). The cause of the two failures is **UNKNOWN**
  (a bad range response, or a damaged scratch file under four parallel fetches, are untested guesses). `curl_resume` now says on stderr how many resumes a transfer
  needed (PopSeq branch `resume-trace-51`, no PR yet) so the next failure can be read.
- The earlier 281 sketches of this set were deleted from the HPCC with the other P7 one-offs; the sketches of the re-fetched 20 are in `bench/T51_retry_p7/` (purgeable).
  The set summary above (255 in analysis, median 38.9 Mb) is **not** redone with the 13 newly recovered runs.
