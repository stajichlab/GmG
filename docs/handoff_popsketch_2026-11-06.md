# Handoff: popsketch, SRA population sets, web page (2026-11-06)

For whoever continues this (a person or another Claude session). Results are in `docs/popsketch.md`; this file holds what is not there: where things are, what is open, what went wrong, and how to read the numbers.

## State at the end of the session
- Merged: #156 (sketch command, bundles, ENA labels, landing/gmg/popsketch pages, SRA scripts), #158 (estimator fix for #157).
- Branch `plan-web-sketches`: `docs/plan_web_sketches.md` (compact hosted sketches), reviewed by Opus, not built.
- Site `gmg.fungalgenomes.org` is deployed (repo variable `PAGES_ENABLED=true` was empty before; the workflow is skipped without it, and a push made before the variable is set is not retried: run `workflow_dispatch`). HTTPS enforcement is still off (owner's setting).
- Open issues: #159 (validate the estimator fix), #155 (real recombining population), #154 (five *B. dendrobatidis* samples with a shoulder). No open PRs.
- Nothing of ours is running on HPCC.

## Where the data is (HPCC)
`/bigdata/stajichlab/jstajich/projects/GmG_popgen/` (work belongs under `~/projects`, not the shared group folder).
- `popsra/<slug>/{runs.tsv,runs_ftp.tsv,sketches/}`: 2,002 sketches (24 GB, 1 Gb of bases each, k=21, scaled 200). Slugs: z_tritici, c_neoformans, c_beticola, p_nodorum, b_cinerea, n_discreta, p_oryzae.
- `popsra/<slug>_out` (default run), `_pcoa`, `_hashpca` (each has `report.html`, `report.pdf`, `figures/`, `pca.tsv`, `tree.nwk`, `clusters.tsv`, `label_check.tsv`, `qc.tsv`, `ani.tsv`), `_sweep` (`--min-floors` 1,2,3,5,8), `_fix` (new binary, estimator check), `<slug>_labels.tsv`.
- `popsra/c_neoformans_bundle.tar.gz` (661 MB, a working web bundle, not hosted).
- `depth/` and `subtest/`: 1, 2, 4 Gb sketches and outputs of the depth tests. The 6 Gb and whole-run sketches were deleted on 2026-11-06 (reproduce with `08_depth_series.sbatch`).
- `popsketch_build/target/release/popsketch` (build), `popsketch_target_fix/` (rebuilt for the estimator test, can go).
- Scripts: `analysis/popgen_sra/01..09*` in the repo.
- Shared hand-off notes with the peer (PopSeq) agent: `/bigdata/stajichlab/shared/projects/Fungi/PopSeq/specifications/GmG_spec_handoff.md`, `GmG_reply_to_004.md`.
- Shared `GmG_runs` (`/bigdata/stajichlab/shared/projects/GmG_runs`): the owner started deleting on 2026-11-06 (about 54 GB left). The published indexes are safe as GitHub release assets (`index-2026.10.3` to `.5`, SHA-256 verified for `.5`). Not published: the finer indexes (`bfd_s1000`, `s2000`, `s50000`, `reps*`), `genomes.list`, `meta.tsv`, cluster tables, `decoys`, `truthset`.

## How to read the results (the traps)
- **Do not use a very deep whole run as the reference.** At 150x and more the old estimator lost the abundance peak (fixed in #158 up to 4 Gb), the cutoff fell to 2, and distances shifted by 0.02. For 15 to 40 Gb *Z. tritici* runs the "whole run" was the wrong answer and 1, 2, 4, 6 Gb agreed with each other (r >= 0.998). Compare depths against each other, and check `coverage_basis` first.
- **Pearson r depends on how diverse the 12 samples are.** The first subsample test (distances mostly < 0.005) gave r 0.62 at 1 Gb for *Z. tritici*; the depth series (distances 0.01 to 0.04) gave 0.998. Mean absolute difference and "mean difference / noise floor" are the safer numbers.
- **"Purity" in `label_check` is misleading** when most labels are singletons or one cluster holds everything (it is 1.0 trivially). Use ARI and the contingency table. Country and study labels are confounded (e.g. 187 of 265 *P. nodorum* are Australian).
- **"No clear cut" is a verdict about the population plus the sketch.** For *Z. tritici* the depth series says 1 Gb resolves distances this size, so it is the population. Not shown for the other species.
- **The one clear cut (*C. neoformans*) is two samples**, SRR34066268 and SRR34066269 (PRJNA1279655). Not identified yet (other species, contaminant, aneuploid?). PC1's 64% is mostly those two. Do not call it lineage structure.
- Samples flagged `unrelated_to_cohort` with a peak at 150x to 240x and a size of 3 to 4 Mb (ERR13824414, SRR26462617) are probably runs dominated by another small genome; unchecked (#159).
- The figures from `_pcoa` and `_hashpca` have not been looked at by eye.

## Operating notes (HPCC and tools)
- Slurm over plain ssh needs `bash -lc`. Interactive runs: `srun -p short` (max 2 h). Do not rebuild a binary while jobs use it; build with a separate `CARGO_TARGET_DIR`. `scontrol update ArrayTaskThrottle` changes the concurrency of a running array.
- `parallel-fastq-dump` (NCBI) gets reads; `sra-stat` 404 means the run is not at NCBI (21 *P. oryzae* ERR runs, skipped by decision, `popsra/p_oryzae_missing.txt`). 7 other failures were transient and succeeded on rerun; concurrency was not shown to be the cause.
- ENA portal API allows CORS (metadata works in the browser); ENA's FASTQ host was too slow and unreliable here; AWS Open Data S3 has no CORS and holds `.sra`. So streaming reads in a browser is parked; sketches must be pre-built.
- Old SRA files are position-sorted, so keep every Nth pair, not the first reads.
- macOS `sed -i` needs `''`; use Python for edits in scripts. `build_site.sh` deletes the output dir: restart any `http.server` afterwards, and new web files must be in its copy list (`ena.js`, `bundle.js` were missed once).

## Claims that were wrong in this session (check the first-hand number before stating)
- Sketch size: said 0.1 to 0.4 MB; measured about 8 MB (whole 3 MB plus two halves of 2 MB). Corrected in the docs and told to the peer agent.
- 28 failed downloads called transient; 21 were permanent 404s.
- Said there was no CNAME and Pages was off; the API showed a CNAME and an approved certificate.
- *C. neoformans* "strong lineage structure": it is a 2-sample outlier cluster.
- Peak counts at 4 and 6 Gb misreported once; corrected in #157.
- Compact sketch size estimate 0.35 to 0.6 MB; the reviewer's arithmetic gives about 0.9 to 1.0 MB (plan v2).
- A threshold lowered to 0.2 produced false peaks; reverted. The accepted rule (0.3 of the hashes seen twice or more in the peak window) is validated only on the samples it was chosen from (#159).

## Decisions made by the owner, not to be revisited by default
- Skip the 21 *P. oryzae* runs. Merging PRs, enabling Pages and HTTPS, DNS: the owner decides (the owner did enable Pages and merge on request this session).
- Flat site layout: `/gmg.html`, `/popsketch.html`, landing at `index.html`.
- ENA is used for metadata only; the web page takes only the user's own files or a hosted bundle.
- Hosting for pre-built sketches: GCP storage or a free S3 on NRP once the sketches exist; plan first, build when HPCC load is lower.

## Next, in order
1. Review the PCA/PCoA figures and trees (cheap); identify the two *C. neoformans* outliers.
2. #159: independent check of the estimator fix (reads of the four changed samples; 6 Gb and whole-run re-test; held-out threshold).
3. Plan `docs/plan_web_sketches.md`: measure the real compact size on one species (gate 1), decide the index-screen handling, then parity on random subsets (gate 2), then HPCC packing when the queue is quiet.
4. Optional: write the results with figures somewhere shareable (the site is public; reports are not published).
5. #155 and #154 when there is a real recombining set and time.
