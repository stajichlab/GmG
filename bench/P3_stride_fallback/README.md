# P3: stride sampling and the windows-with-fallback path, tested on the biased runs (HPCC)

The phase 1 pilot (`bench/P2_phase1_pilot/`) found that contiguous windows are biased on position-sorted runs: two of ten
*C. auris* runs came out as a 2.2 and a 2.3 Mb genome. PopSeq's `popseq fetch` now samples the coverage tier by evenly
spaced reads by default (`--sampling stride`), keeps `windows` as an opt-in, and falls back from windows to stride when the
hook's check fails (hook exit 75). `popsketch check` is that check. This folder is the first test of both on real data:
the two biased runs (ERR2300774, ERR12543635) and one sound control (DRR129824), 25x, `--require-pairs`, the sketch hook,
two arrays of three tasks at once (declared concurrency 6). `summarize.py` builds `summary.md` from the raw files here.

## Result (n = 3 runs of one species; the same three runs in each array)

See `summary.md`. In short:

| | Stride (default) | Windows with the check and fallback |
|---|---|---|
| ERR2300774 (biased) | ok, 11.9 Mb, 19x | windows failed the check (2.2 Mb, 106x, 4.2 times the expected coverage); resampled by stride; ok, 11.9 Mb, 19x |
| ERR12543635 (biased) | ok, 12.2 Mb, 19x | windows failed the check (2.3 Mb, 107x, 4.3 times); resampled by stride; ok, 12.2 Mb, 19x |
| DRR129824 (sound) | ok, 11.9 Mb, 20x | windows passed the check (11.9 Mb, 20x); no fallback |

- **The biased runs now give the right genome.** The estimated size went from 2.2 and 2.3 Mb to 11.9 and 12.2 Mb (genome
  12.25 Mb), and the k-mer coverage from 106 and 107x to 19x. Coverage reached 25.000x for both (24.985x for DRR129824).
- **The check separated the cases cleanly:** size ratio 0.97 to 0.99 and coverage ratio 0.76 to 0.80 for every sound sample,
  0.18 and 0.19, and 4.2 and 4.3 for the two biased ones. The thresholds (0.5 and 2.5) fall in the gap, but the gap is
  n = 5 samples of one species.
- **The fallback did not download again.** The log has, for each biased run, one `ok` row with the transfer (283 MB and
  671 MB) and a second `ok` row with no bytes whose `error` says "resampled evenly spaced from the same transfer"
  (`downloads_p3.csv`). The windows attempt, the check, the resample (65 s and 62 s) and the second check all
  happened inside one task.
- **No reads were kept** (the hook sketched and the reads were deleted).

## What it costs

| | Stride | Windows (phase 1) |
|---|---|---|
| Result seconds, the same three runs | 646, 313, 360 | 166, 90, 89 |
| Ratio | 3.9, 3.5, 4.0 times | |
| Peak memory (task) | 1.9 to 3.1 GB | 0.6 to 2.2 GB (phase 1, all ten) |

Stride converts the whole run: 3.5 to 4.0 times the seconds of windows here. A sorted run taken by windows with the fallback
costs about what stride costs (341 and 354 s against 313 and 360 s), because the windows attempt is cheap; a run that is
not sorted costs the windows time (161 s). Whether windows with the fallback is cheaper than stride therefore depends on the
share of sorted runs, which is unknown: 2 of 10 in the one species measured.

## Limits

- **n = 3 runs, one species, one genome size (12.25 Mb).** The thresholds of `popsketch check` are not tested on any
  other species, on larger or polyploid genomes, or on contaminated samples (a contaminated sample could fail the size
  or coverage rule for another reason and be resampled needlessly).
- **Timings are noisy.** Six tasks ran at once and may have shared nodes (PopSeq's concurrency sweep found that the
  speed per fetch depends on how many fetches share a node), and the phase 1 times were measured at concurrency 5.
- **That the two runs are position-sorted is still inferred** from the match between the sketch and `25 / depth x genome`;
  the read order of the `.sra` was not inspected.
- Peak memory for stride was 3.14 GB at most (n = 3); the 8 GB requested was more than needed, 4 GB is tight, 6 GB is safe.

Files: `manifest_test.csv`, `stride/` and `windows/` (`results.tsv`, `sampling_check.tsv`, `hook_times.tsv`), `sacct.tsv`,
`downloads_p3.csv`, `summarize.py`, `summary.md`. On HPCC the work is in
`~/projects/GmG/bench/P2_phase1_build/bench/P2_phase1_pilot/work_stride` and `work_windows`.
