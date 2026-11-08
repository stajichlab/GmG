# Data problems

Known problems in input data, and how each one is handled. One row per problem; add new ones at the end. Runs that are left out of an analysis are in `analysis/popgen_sra/excluded_runs.tsv`, which `04_run_species.sh` reads (the sketch stays on disk; it is just not passed to `popsketch run`).

| Run(s) | Set | Problem | Handling | Status |
|---|---|---|---|---|
| SRR34066268, SRR34066269 | *C. neoformans* (SRA population set), both PRJNA1279655 | Outlier pair. They make the only clear cut (262 vs 2, gap/floor 26.9) and most of PC1 (64%), so they hide any structure in the other 262. What they are is not known: a different species, a contaminant, or an aneuploid isolate are all possible. | Suppressed: listed in `excluded_runs.tsv`; rerun the species without them. | Suppressed; cause not investigated. To check: k-mer containment against *C. neoformans* H99 and a few *Cryptococcus* species, and the abundance histogram. |
| 21 *P. oryzae* ERR runs (ids in `popsra/p_oryzae_missing.txt` on HPCC) | *P. oryzae* (SRA population set) | Not downloadable: `sra-stat` returns 404 at NCBI, and the runs were not found at ENA either. Sequence cannot be found. | Dropped for good (2026-11-07). The set has 251 sketches, not 272. The page catalog should list them as excluded with this reason (see `docs/plan_web_sketches.md`). | Closed |
| SRR27558955 | *C. neoformans* (SRA population set), PRJNA1065146, labelled India: Varanasi, strain NK1 | Second genome in the reads: the genome-size estimate is 45 to 120 Mb at every depth from 1 Gb up, and the abundance peak is not the target's. | Not suppressed yet; the estimator in #158 does not handle mixed-genome reads. Flagged by the size and `unrelated_to_cohort` checks where it applies. | Open, not suppressed (you asked only for the two above); say if it should go in `excluded_runs.tsv` too (#159) |

Related, not a data problem: SRR10298347 estimates a wrong size at 4 Gb only (estimator, #159).
