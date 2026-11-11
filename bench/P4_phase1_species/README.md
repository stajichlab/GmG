# P4: phase 1 of the PopSeq panel plan, one full species (*Brettanomyces bruxellensis*, taxid 5007)

All 300 panel isolates of one species, fetched at 25x by stride and sketched, a few fetches at a time per node (plan:
`docs/plan_popsra_scaleup.md`, phase 1). Status: **prepared, not launched**.

## Why this species
Cheapest of the 46 species with 300 isolates and a genome size that I checked, apart from *K. lactis* (172 of its 300 isolates
are under 25x): 344 Gb of runs in all, a genome of 13.2 Mb, 1 isolate under 25x and 3 over 200x, 14 studies. *C. auris*, the
species used in the 10-isolate pilots, is 1,071 Gb (3 times larger; 143 isolates over 200x).

## What is in the manifests
- `manifest.csv`: 299 runs, 296 Gb. `manifest_deferred.csv`: SRR3927079, 47.9 Gb (3,641x). Stride converts the whole run at
  about 0.19 to 0.40 s per Mb (n = 3 runs), so that run alone would take 2.5 to 5.3 hours and could exceed the job's 4 h limit
  and take its chunk down with it. It is deferred, to be taken by `--sampling windows` with the check, or by a longer job. The
  next largest run is 7.7 Gb (0.4 to 0.9 h).
- Built with `bench/P2_phase1_pilot/select_phase1.py <panel metadata> 5007 300 manifest.csv 0 1e12`, then split at 10 Gb.

## Expected cost (arithmetic from the phase 1 pilot, not measurements at this scale)
- 296 Gb x 0.19 to 0.40 s/Mb = 16 to 33 task-hours of stride conversion, and 83 to 192 GB to transfer
  (0.28 to 0.65 bytes per base, n = 3). With 20 fetches at once, about 1 to 2 hours plus queue time.
- 10 jobs (chunks of 32 at 4 fetches each), at most 5 at once: 5 x 4 = 20 fetches at a time, 4 per node.

## How to run
See the header of `scripts/popsra_fetch_chunks.slurm`; failures are turned into a retry manifest by
`scripts/popsra_retry_list.py MANIFEST $WORK/results.tsv -o retry.csv`. Then `popsketch run --from-sketches $WORK/sketches
--collapse mean --collapse-within study --labels ...` for the species result.

## Gate (from the plan)
At least 95% of runs `ok`, and the per-run cost known (wall time, memory at 4 fetches per job, the share routed to ENA, the
share of runs the sampling check fails, the transfer rate per fetch at 4 per node).
