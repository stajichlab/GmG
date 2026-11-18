# S52: the species check in the fetch hook, end to end (2026-11-19)

Two real fetches (`popseq fetch`, tier coverage 25, sampling auto, PopSeq 0.2.0, job 32174455 on stajichlab, 8 CPUs, 2 at a time) through `scripts/popsra_sketch_hook.sh`
with the real `popsketch species-match` and `bfd_s10000.gmgi` (`SPECIES_MATCH=refuse`). Files: `results/` (`exits.tsv`, `species_match.tsv`, its `.meta.json`).
`run.sh` + `manifest.csv` re-run it; everything it writes is under `out/` (purgeable).

| run | filed | exit | verdict | evidence | sketch |
|---|---|---|---|---|---|
| SRR30534096 | *C. neoformans* | **14** (`hook_failed`, hook exit 66) | `mismatch`, `exclude_default` yes | *Aspergillus*, best hit 42% of the hashes, 0.2% in the filed genus | moved to `sketches/refused/` |
| DRR829882 | *C. neoformans* | 0 | `match`, `exclude_default` no | *C. neoformans*, 96% in the filed genus | stays in `sketches/` |

The `species_match.tsv` has one header and one row per run. The known mislabel is refused and the clean control is kept, with the real binary and the small index inside a real fetch job.

What this does not show: the unchecked paths (no index, no filed name) and the six-at-once race are covered by `scripts/test_popsra_sketch_hook.sh` with a stubbed `species-match`, not here; two runs of one species, one index, one node.
The job log has a cosmetic line from existing code (`popsra_fetch_chunks.slurm` line 83: `exits.tsv: No such file or directory`, from reading the line count of a file that does not exist yet); it is not from this change and was not fixed.
