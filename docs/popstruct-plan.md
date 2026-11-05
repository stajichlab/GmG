# Plan: `popsketch` (audit a set of read files; coarse population structure as a secondary output)

Status: decisions agreed 2026-11-04 after a question-by-question review of the first draft. Issue #15. Written after the
Coccidioides and *B. dendrobatidis* explorations in `analysis/popgen_cocci/` (scripts 00-17). Not implemented.

## Goal and headline

Give the tool a set of FASTQ (or FASTA) files and, optionally, a label file. It returns:

1. **The audit (headline).** Non-target organisms, contaminated, half-content and low-coverage samples, and label
   conflicts, each with a stated reason and the numbers behind it. This is what the explorations actually found: 3
   non-*Coccidioides* read sets, 2 label conflicts and 26 half-content *B. dendrobatidis* samples. (An earlier claim of 31
   contaminated Coccidioides samples was mostly an abundance-cutoff artifact; see Findings.)
2. **Coarse population structure (secondary).** Hash PCA, clusters and a first-pass tree, always shown together with the
   resolution limits of the data. The report opens with the audit; the structure follows.

Success means "it finds the known problems in the two real datasets", not "the tree looks right".

## Decisions

| # | Decision | Basis |
|---|---|---|
| 1 | Audit is the headline; PCA, clusters and tree are secondary and carry the resolution advisor. | Sketches matched SNP structure on the top axes only; distance agreement among close pairs was 0.05-0.33. |
| 2 | QC flags are **cohort-relative** (robust spread around the cohort median; a minimum floor so a tight cohort is not over-flagged). The index is optional extra evidence, never required. If more than about 30% of samples are flagged the tool warns that the cohort median cannot be trusted. | Works for any organism, including ones the index covers poorly. Every flag prints its reason and numbers. Reconfirmed 2026-11-04 after the abundance-spectrum experiment: the index stays optional, the cohort-only path reports heavy contamination only, and the report states that lighter contamination (1-7%) is not assessable without `--index`. |
| 3 | Flagged samples are **excluded by default**. `--flagged exclude\|mark\|decontaminate` selects the mode. Decontamination is opt-in, drawn with its own marker in every plot, and kept out of the resolution advisor. | Excluding fixed the Bd PCA (cosines 0.999 vs SNP PCA). The decontamination vote removes private sequence and pulls samples toward the cohort, so its size result is not an accuracy result, and the size gain it showed was largely an abundance-cutoff artifact (see Findings). |
| 4 | Label file has **several columns**; the first is primary. A blank or `NA` cell means unlabelled. The conflict check needs 3 or more members in a group. Each column gets its own colors, ARI and conflict check. | Both datasets needed more than one label (species/lineage, plus source). Labels enter after the distances, so extra columns are cheap. |
| 5 | **R1 and R2 are merged** per sample with the same grouping rules as `gmg run`; `--r1-only` is a fast mode. | Consistent with `run`; extra coverage reduces false low-coverage flags. |
| 6 | Input cap is **fixed bases per sample** (R1 and R2 combined), default about 1 Gb, documented as a default to adjust. The audit reports achieved coverage and warns when its spread across the cohort is large. | Equal reads did not give equal coverage (2x to 59x with a 6M-read cap). Target-coverage stopping is deferred. |
| 7 | Design and test for about **1,500 samples** with exact all-pairs; refuse above a configurable `--max-samples`, suggesting scaled 1000, a subsample or a split. | About 5e10 comparisons at n = 559 took minutes; n = 1,500 stays tractable; an inverted index is not needed yet. |
| 8 | **Hash-presence PCA is the default ordination**; PCoA of 1-ANI stays as `--ordination pcoa`, written to a file but not in the headline figure. | Hash PCA tracked SNP PCA better on all four datasets (table below). |
| 9 | **Clusters:** cut UPGMA at the largest merge-height jump only when that jump is at least 3 noise floors (`--min-floors`); otherwise report "no single clear cut" with the top 3 candidate cuts and their jumps. Everything merging below the floor is one population (k = 1). `--k` and `--cut` override. The gap ratio gap1/gap2 is reported with its probability under a structureless coalescent population, but is not the gate. **Changed from the first agreement (ratio of at least 2) on new evidence:** that ratio is at least 2 in about 62% of structureless populations (25% at 5, 10% at 10; `analysis/popgen_cocci/26_gap_null.py`), so it cannot tell structure from the natural asymmetry of a coalescent tree. | The jump must exceed the measured technical noise to count. On the real sets: both Cocci species 6.2 floors (clear, k=2), Bd 0.9 floors (no clear cut; k=5 is the top candidate, matching the five coarse groups found by hand). The threshold of 3 floors is a judgement from these two sets. |
| 10 | **Tree:** NJ on mash distance with unsupported branches collapsed to polytomies; unrooted by default, `--root ID` available. The collapse floor is the split-half noise floor (see below). | Within-lineage trees were combs; a fully resolved tree would imply resolution the data lacks. Midpoint rooting failed whenever an outlier was present. |
| 11 | **Resolution limit** is a split-half noise floor measured per run: alternate reads of 20 samples go into two sketches and the floor is the 90th percentile of the distance (1 - ANI) between the halves. It collapses NJ branches below it and drives the resolution advisor and the gap rule. With too few read sets to split, a fixed floor of **0.001** is used. **Resolved, with a change:** the fallback is 0.001, not "about 10 differing hashes" (10 hashes is about 3e-6 in distance, three orders of magnitude below the measured floors of 0.0002 to 0.002). | Experiment 2026-11-04: the 90th-percentile floor predicts that Bd lineages 2/3/4 are indistinguishable (99% of those pairs below it) where the median does not (25%). The check against SNP distance was inconclusive; the floor is pessimistic (each half has half the coverage). |
| 12 | **Home:** a new crate `crates/popsketch` (lib plus bin), added to workspace `members`, after a prep PR that extracts the shared input handling into a small library crate. | The code to reuse (`sample_key`, `group_reads`, `read_sheet`, discovery, `--resume`) is private inside the `gmg` binary crate (`gmg-cli/src/run.rs`) and cannot be imported. Copying it would drift. |
| 13 | **Tests:** a deterministic generator written in Rust inside the tests plants a population (about 60 samples, 3 groups, known tree) plus problems: 2 low-coverage, 2 contaminated (10-50% foreign reads), 2 half-content, 1 foreign genome, 2 swapped labels. Pass: every planted problem is flagged, no clean sample is, swaps appear as conflicts, planted groups come back. The Cocci and Bd runs on the HPCC remain the check on thresholds. | Real data cannot ship in a public repo. A fixture proves the logic, not that the thresholds are right. |

## v1 scope

In v1:
- shared input handling, merged R1 and R2, the per-sample base cap;
- cohort-relative QC flags with reasons and the three flag modes;
- the optional `--index` organism screen;
- containment ANI and mash distances, exact all-pairs, up to about 1,500 samples;
- hash PCA, with PCoA as an option;
- largest-gap clusters and the collapsed NJ tree;
- the multi-column label file with conflicts, ARI and a contingency table;
- the resolution advisor;
- an HTML report plus TSV and JSON outputs;
- the simulated test fixture.

After v1: external comparison (`--compare-tree`, `--compare-pcs`: exact-clade test, principal angles, Procrustes, distance
Spearman); DAPC; Markdown and PDF reports; a browser (wasm) build; target-coverage stopping; an inverted index beyond 1,500
samples.

## What exists and can be reused

- `gmg-core`: FracMinHash sketching with abundances, `.gsk` files, `auto_min_abund` (`coverage.rs`), containment and Jaccard
  helpers, classification. Only dependency is a hash crate, so it stays wasm-friendly and must stay that way.
- `gmg-cli/src/run.rs` (binary crate, not importable): input discovery, R1/R2/lane grouping, sample sheets, `--resume`
  (keyed on path and size), the report writer. Reuse requires the prep PR (decision 12).
- `gmg-report`: hand-written SVG charts (ANI histogram, heatmap with a fixed 34 px cell, bars). A scatter plot with hover,
  tree drawing and a scalable heatmap are new work.
- Not in the workspace: linear algebra, tree building, clustering. These go in the new crate; any heavy linear-algebra
  dependency sits behind an optional feature, never in `gmg-core`. `rayon` lives in `gmg-cli`, not `gmg-core`.

## Target interface (flags are proposals; names may change)

```
popsketch run reads/ [--sample-sheet S.tsv] [--labels labels.tsv] -o out/
             [--scaled 200] [--max-bases 1000000000] [--r1-only] [--min-cov 4] [--max-samples 1500]
             [--flagged exclude|mark|decontaminate] [--index bfd.gmgi]
             [--ordination hash-pca|pcoa] [--k N | --cut D] [--gap-ratio 2] [--root ID] [--tree-floor auto|N]
```

Outputs: `qc.tsv` (flags, reasons, achieved coverage, abundance cutoff, estimated size), `ani.tsv`, `mash.tsv`, `pca.tsv`,
`clusters.tsv`, `tree.nwk`, `label_check.tsv`, `popsketch.json`, `report.html`. Resume covers sketching only; the abundance
cutoff, QC, flags and everything after depend on the whole cohort and are recomputed.

## Stages

1. **Sketch.** One sketch per sample (R1 and R2 merged), scaled 200 (1000 derivable), abundances tracked, capped at the
   per-sample base limit. Two sketches per sample (alternate reads) are kept for the split-half noise floor. The abundance
   cutoff is `gmg-core`'s `auto_min_abund` (the valley between the error tail and the coverage peak; 2 when there is no
   peak). It is a QC output and is exposed. A rule that returned 2 for almost every sample (the Python prototype's) inflated
   sizes and distances and must not be used.
2. **QC.** Flag by cohort-relative size, coverage and shared-hash fraction; low coverage (below 4x; 14 Cocci and 1 Bd
   samples); half-content samples (26 Bd at 11.8 Mb against about 21, robust to the cutoff choice); samples unrelated to the
   cohort (NM_7898, CA10 and CA12 shared almost no hashes). With `--index`, name the organism (CA10 and CA12 were
   Lecanicillium, NM_7898 Histoplasma). **Contamination (cohort-only path): heavy cases only.** A cohort-relative size test
   catches heavy contamination at low coverage; an abundance-spectrum "shoulder" test catches a contaminant whose own coverage is
   several-fold; neither sees light contamination (see "Contamination: abundance-spectrum experiment"). The report must state
   that detection limit. Light contamination (1-7%) is found only with `--index`, via the existing `mixture_suspected` /
   `second_taxon`. Decision 2 (the index is never required) stands, with that stated limit.
3. **Distances.** Containment ANI (C^(1/k), C = shared / min size) and mash distance from Jaccard, sorted-hash merges with
   rayon. ANI is labelled an estimate.
4. **Hash PCA.** Presence/absence over hashes seen in at least 2 samples and not in all. Because a private or error k-mer adds
   a variance direction only to its own sample, the Gram-matrix route needs a second pass restricted to this frequency-filtered
   set; it cannot simply reuse the all-hash intersection counts used for distances. Report each PC's correlation with
   estimated size and coverage and warn when it is high (Bd k-mer PC1 was a size artifact).
5. **Clusters and tree.** Decisions 9 and 10.
6. **Labels.** Per column: colors, ARI, contingency table, purity; per-sample conflicts from median ANI to each label plus
   nearest-neighbor label (this found the two Coccidioides swaps).
7. **Resolution advisor.** From the noise floor: the fraction of sample pairs closer than the floor, and a statement that
   structure below it is not supported, with a recommendation to use a SNP method. Decontaminated samples are excluded from it.

## Noise floor (split-half): experiment result (2026-11-04)

While sketching, alternate reads go to two sketches of the same sample. The distance (1 - containment ANI) between a sample's
two halves is the technical noise. Run on 12 Cocci and 12 Bd samples (halves of up to 4M reads, gmg-core cutoff, scaled 200;
`analysis/popgen_cocci/20_splithalf_*.sh`, `21_splithalf_analyze.py`). Halves have less coverage than a full sketch, so the
floor is pessimistic; no rescaling has been derived.

| | split-half distance, median | 90th percentile | max |
|---|---|---|---|
| Cocci (12) | 0.00018 | 0.00204 | 0.00326 |
| Bd (12) | 0.00047 | 0.00163 | 0.00392 |

- **Prediction 1 holds with a high quantile, not the median.** Bd lineages 2/3/4 (within-group median distance 0.0007): 99%
  of pairs fall below the 90th-percentile floor but only 25% below the median floor. Pairs between {2,3,4} and the other
  groups (median 0.0049) are essentially never below either floor. Cocci within-species pairs (median 0.003-0.005) are above the
  floor (0.2-3% below).
- **Prediction 2 is inconclusive.** Agreement with SNP-PC distance among pairs below vs above the floor: Bd median floor 0.39
  vs 0.83 (supports), Bd 90th-percentile floor 0.77 vs 0.34 (reversed), Cocci 0.20 (n=80) vs 0.83 at the median floor and 0.82
  vs 0.83 at the 90th percentile. These correlations are confounded by group structure (above-floor pairs fall into a few
  discrete group distances), so they do not separate the floor from that effect.
- The floor is heavy-tailed: a few low-coverage halves have floors 10x the median, so the choice of statistic matters.
- Twelve samples per dataset; one cutoff; a SNP-PC proxy for truth.

Consequence for decision 11: the split-half floor is usable if defined conservatively (high quantile, or per pair from the two
samples' own floors), and it could be estimated from a random subset of about 20 samples instead of all of them. Whether it
should replace the fixed fallback needs a decision, and a better check than SNP-PC correlation (for example simulated data
with known branch lengths from the test fixture).

## Contamination: abundance-spectrum experiment (2026-11-04)

Question: can the k-mer abundance spectrum alone, with no index, flag mixed samples? A contaminant at read fraction f has its own
coverage of about c * f / (1 - f), where c is the host coverage; at 2.6% and 43x that is about 1x, indistinguishable from errors.
Two cohort-relative statistics (robust z within coverage bins; `analysis/popgen_cocci/22_abundance_spectrum.py`): the
*shoulder* (hashes at abundance 0.1c-0.6c as a fraction of the genome mass) and its *excess* over a Poisson(c) expectation.
Samples with c below 6 are not assessable by the spectrum. Truth: `gmg run --index` mixture calls (`second_explained_fraction`),
which are another method's estimates, not ground truth. Existing results: `~/projects/TE_fungi_copy/results/{cocci,bd}_gmg/out`.

| | Cocci (558) | Bd (334) |
|---|---|---|
| index-confirmed mixtures | 20 (7 at 5% or more) | 1 (6.3%) |
| assessable by the spectrum (c >= 6) | 15 of the 20 | 1 of 1 |
| detected by the spectrum (z > 3) | 2: B12526 (20%, c=41, z=11.7), B13956 (1%, z=5.3, also size z=7.5) | Sample_55930 (z=3.6) |
| mixtures at 1-3% | 0 of 13 detected (AUC 0.50, chance) | |
| not assessable (c < 6) | 5: TX7, M192, M216, M158, NM_0071 (31-46%), all caught by the size test (z 88-100) | |
| M223 (7%, c=16, own coverage about 1.2x) | missed (z=2.2) | |
| false positives (z > 5 among non-mixtures) | 7 of 492 (1.4%); at z > 3, 11 of 492 | 7 of 324 (2%); at z > 3, 19 of 324 |

- Size and spectrum are complementary: size catches heavy contamination at low coverage (5 of 7 mixtures at 5% or more); the
  spectrum catches a contaminant at several-fold coverage (B12526, which size misses at z=0.9). Together 6 of 7 at 5% or more
  and 0 of 13 at 1-3%.
- The 7 Cocci "false positives" are not all false: two are cutoff misfires (COCPO_103716, NM_4297, z in the hundreds, 0.1-0.2 Mb
  of hashes) and the rest are unexplained (B16536, NM_3006, VFC090, NM_4708, B11199); the index may simply not see a contaminant
  that is not in the reference or is split across taxa. In Bd the z > 5 list includes the half-content samples. Not resolved.
- Statistical power is poor: only 2 assessable Cocci positives at 5% or more, and one Bd positive. This sets a detection limit;
  it cannot set thresholds.

Consequence: the cohort-only path reports heavy contamination (roughly 30% or more at low coverage, or about 10% or more at high
coverage) and says plainly that lighter contamination is not assessable; `--index` adds the light cases.

## Thresholds that remain judgement calls

The MAD multiplier (5) and size floor (3% of the median) for QC flags; 3 floors for a clear cut; the 1 Gb cap; the 90th percentile as
the floor statistic and 0.001 as its fallback; the shoulder z of 5. They were set from the Coccidioides and Bd sets (and the
simulated fixture), not tuned on independent data, and are all options or documented constants.

## Build order and status

0. Experiments (done 2026-11-04). The existing abundance cutoff works except for four sample histograms; the guard is in
   PR #146. The split-half floor is usable as a high quantile.
1. Prep PR: `gmg-input` crate extracted from `gmg run`, behavior-preserving (PR #147).
2. `popsketch` crate (this PR): sketching with a per-sample base cap and split-half sketches, cohort QC with reasons, the optional
   index screen, exact distances, hash PCA (Gram matrix with its own eigensolver) and PCoA, UPGMA with the noise-floor cut,
   collapsed NJ tree with Newick output and rooting, multi-column labels with conflicts and ARI, the resolution advisor, the
   HTML report, `--flagged exclude|mark|decontaminate`, `--max-samples`. 48 unit tests and 3 end-to-end tests on simulated reads
   with planted problems (a deterministic Rust generator): exactly the planted problems are flagged, the planted groups and the
   swapped label are recovered, a clonal set is one population.
3. Validation on the real sets (`analysis/popgen_cocci/27-29`; numbers in `docs/popsketch.md`): done for Cocci (559) and Bd (334).
4. Not done (deferred, as agreed): external comparison (`--compare-tree`, `--compare-pcs`), DAPC, Markdown and PDF reports, a
   browser build, target-coverage stopping, an inverted index beyond 1,500 samples; a tutorial.

## Findings the design rests on

**Correction (2026-11-04).** The first explorations used a Python abundance cutoff that returned 2 for almost every sample.
Rerunning with gmg-core's `auto_min_abund` (no decontamination) changed several results, listed here. Numbers below are from
the reruns unless marked.

- Coccidioides: species gap in ANI (within about 0.996, between about 0.980); both species and the 19-sample immitis cluster
  are exact SNP clades (that test used the earlier cluster assignments); posadasii k=3 supported for one split only (12 of 262
  samples discordant); two `samples.csv` label problems (not re-checked with the corrected cutoff: the rerun applied the
  earlier relabelling).
- **Retracted:** "31 contaminated Coccidioides samples" and the decontamination gain (32.6 to 27.4 Mb). With the proper cutoff
  the cohort is tight (25.9 Mb median, 5-95% 25.5-26.7 at 8x or more); 25 of the 33 flagged samples are no longer outliers. **This
  retraction overreached for the rest:** five of the 8 samples that remain (TX7, M192, NM_0071, M158, M216, at 2.2-5x
  coverage and 53-57 Mb) are index-confirmed mixtures (31-46% of hashes explained by a second taxon, mostly Histoplasma), so
  their size inflation was real contamination. Of the 20 index-confirmed mixtures in the Cocci set, the nine strains first
  reported are at normal size (25.4-28.4 Mb), so size misses light contamination (see the experiment below).
- Bd: nearly clonal, all pairs ANI at least 0.98; 5 coarse groups recovered (merged after seeing the ANI table, so not
  independent); the 26 half-content samples stay under 17 Mb with either cutoff; cleaned k-mer PCA matches the genome-wide VCF
  PCA (top-3 cosines 0.998, 0.997, 0.934; distance Spearman 0.968).
- Distance Spearman against SNP-tree distance among close pairs was about 0.05-0.33 in the earlier runs (not repeated with the
  corrected cutoff): sketches give coarse structure only.
- Third-axis agreement is unstable: for immitis and posadasii the third cosine fell from 0.74/0.69 to 0.03/0.14 with the
  corrected cutoff. Quote the top two axes only.

**Hash PCA vs PCoA of 1-ANI, each against SNP PCA (same shared samples; `analysis/popgen_cocci/16_pcoa_vs_snp.py` with
`SUF=_rc`, corrected cutoff):**

| Dataset | cosines top 3, hash PCA | cosines top 3, PCoA | distance Spearman, hash PCA | distance Spearman, PCoA |
|---|---|---|---|---|
| Cocci, both species (n=423) | 0.996, 0.986, 0.091 | 0.987, 0.096, 0.019 | 0.935 | 0.847 |
| *C. immitis* (n=147) | 0.967, 0.884, 0.028 | 0.822, 0.227, 0.155 | 0.870 | 0.521 |
| *C. posadasii* (n=262) | 0.979, 0.963, 0.135 | 0.958, 0.724, 0.044 | 0.901 | 0.791 |
| Bd cleaned (n=256) | 0.998, 0.997, 0.934 | 0.991, 0.965, 0.073 | 0.968 | 0.825 |

Hash PCA is better than PCoA on distance agreement in all four, with the earlier cutoff as well (decision 8 stands).

**Largest-gap cluster rule on UPGMA merge heights (`analysis/popgen_cocci/17_gap_rule.py`, `SUF=_rc`):**

| Dataset | rule picks | gap1/gap2 | comment |
|---|---|---|---|
| Cocci, both species (n=541) | k=2 | 8.9 | the two species |
| *C. immitis* (n=171) | k=2 | 69.2 | consistent with the SNP cluster |
| *C. posadasii* (n=370) | k=2 | 6.0 | consistent with one deep SNP split |
| Bd cleaned (n=307) | k=5 | 1.96 | matches the five coarse groups (defined after the fact); just under the threshold of 2 |

With the earlier cutoff the ratios were 3.0, 3.3, 15.8 and 1.1 (Bd picked k=6).

## Caveats

- Everything before 2026-11-04 afternoon used a faulty abundance cutoff; only the reruns marked as such use the correct one,
  and the SNP-tree clade tests (script 08, immitis cluster and posadasii k=3) were not rerun. The PCA comparison (script 16)
  was.
- All results come from two datasets, one clonal and one with clear species structure. The method is untested on recombining
  populations and on mixtures.
- Several choices were made after seeing the data: the 17 Mb size cut for Bd, the merged Bd groups, k for k-means, the gap
  ratio of 2 (four datasets only).
- The SNP comparisons are not independent of the sketches: the SNPs came from mapped reads, likely the same read sets.
- Missing evidence: a negative control (shuffled labels, a single-population null), sensitivity to scaled, cutoff and sample
  size, variance over repeated read subsamples. Cosines of 0.98-0.999 are for the top few axes and may mostly reflect the
  species or clade split; quote the number of axes and variance share with them.
- The two label swaps are supported by k-mer neighbours and, for 485B-1, by the SNP tree; B3476 is not in the SNP set, and
  neither is confirmed against the sample metadata.

## Review notes (2026-11-04, Fable review; claims checked against the code)

Folded into the decisions above: separate solvers for PCA and PCoA with linear algebra behind an optional feature in the new
crate; adding the crate to workspace `members` (the root `Cargo.toml` omits `gmg-py`); scale target; resume covers sketching
only; smaller report-engine reuse than first stated; the PC-versus-size/coverage check; the gap-ratio null for clusters;
marking decontaminated samples; exposing the abundance cutoff; moving the simulated fixture earlier; putting the CLI before
the report charts; deferring external comparison.

Checked and qualified:
- The reviewer said the Rust cutoff has a known failure on noisy high-coverage histograms and that the fix is to port the Python
  rule. Tested 2026-11-04: not borne out. The Rust rule works on 556 of 559 Cocci and 333 of 334 Bd sketches (it misfires on 4,
  needing a guard); the Python rule is the one that is wrong (it returns 2 for almost every sample). The reviewer was right that the two rules differ.
- The suggested reuse of `classify.rs:12` (`READS_SPECIES_MIN_ANI`) does not apply: it lowers a species cutoff for
  classification and does not correct an ANI estimate.
- That min-size normalisation inflates the between-species ANI of 0.978 is a hypothesis to test, not a finding; the two
  species have similar genome sizes.
