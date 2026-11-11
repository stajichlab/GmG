# P1: laptop pilot of `popseq fetch --tier coverage`

Three species from the PopSeq panel `pop_sra_ge10_cap300`, one run each, 25x, ENA only (the laptop
has no SRA toolkit), three tasks at once (declared concurrency 3). Results and what they show:
`docs/plan_popsra_scaleup.md`, section "The laptop pilot".

| File | What |
|---|---|
| `select_pilot.py` | picks the runs from the panel and writes `plan.csv` (ENA values only: with no NCBI counts the SRA routes fail at once and ENA does the work) |
| `plan.csv` | the plan: six candidate runs, three were fetched |
| `run_pilot.sh` | the three fetches in parallel, from `$SCRATCH` |
| `results.tsv` | one line per fetch attempt (2 ok, 1 verify_failed, 1 retry stopped as transient) |

Not kept (untracked, purgeable): `popseq/` (a clone of PopSeq at e6c5ba0), `reads/` (the two
verified outputs, about 400 MB), `work/` (logs and scratch). To purge: `rm -rf work reads popseq`.
Run with `uv` from this folder after `git clone https://github.com/stajichlab/PopSeq.git popseq`.
Data processing was done on a laptop here by the owner's instruction for this pilot, an exception
to the rule that analysis runs on HPCC (HPCC was unreachable at the time).
