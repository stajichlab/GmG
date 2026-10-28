# GmG work log

Chronological record of what was done, why, and what was decided. Newest at the bottom. Overnight run authorised 2026-10-27 (ssh, sbatch, srun, cargo, git push preapproved; self-merge allowed; HPCC run dir `/bigdata/stajichlab/shared/projects/GmG_runs/`; partitions short/preempt/exfab OK).

## 2026-10-27 (evening)
- Repo created (`stajichlab/GmG`, private), initial scaffold committed (`cae9708`). `gmg-core` FracMinHash is hash-identical to the `sourmash` crate on A1163 and Af293.
- HPCC has cargo 1.96 (x86_64), seqtk module, no read simulator installed.
- Overnight plan: stages 1-7 (see issues). Work on branch `overnight`, merge to `main` at stage gates.
- Stage 1-2 code done (commit on `overnight`): `gmg-core` gained `skfile` (.gsk sketch files, with abundances), `index` (inverted hash->genome index built in hash-range buckets; save/load; `search`), `call` (top-hit/vote logic, hierarchical per-rank calls; blank ranks are skipped, not failed), `downsample`, abundance sketching. `gmg-cli` has `sketch`, `build`, `guess`, `info`, `bench` (loo/los/log, sharded for SLURM). 10 unit tests + sourmash oracle pass.
- Decision: one sketching pass at scaled=1000, coarser indexes (2000/5000/10000) are built by downsampling the sketches (`gmg build --scaled`), so no re-sketching is needed to compare scaled values.
- Decision: calibration of per-rank thresholds is done offline from `gmg bench` TSVs (features: top shared, ANI, vote consistency); the defaults in `Thresholds::default()` are placeholders.
- Smoke test on the 15 pilot genomes: sketch 15 genomes 1.2 s wall (5.5 s CPU, 0.37 s/genome), 2.3 MB index, self-guess ANI 1.0, LOO nearest neighbour correct.
