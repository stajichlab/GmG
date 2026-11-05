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
| 2 | QC flags are **cohort-relative** (robust spread around the cohort median; a minimum floor so a tight cohort is not over-flagged). The index is optional extra evidence, never required. If more than about 30% of samples are flagged the tool warns that the cohort median cannot be trusted. | Works for any organism, including ones the index covers poorly. Every flag prints its reason and numbers. |
| 3 | Flagged samples are **excluded by default**. `--flagged exclude\|mark\|decontaminate` selects the mode. Decontamination is opt-in, drawn with its own marker in every plot, and kept out of the resolution advisor. | Excluding fixed the Bd PCA (cosines 0.999 vs SNP PCA). The decontamination vote removes private sequence and pulls samples toward the cohort, so its size result is not an accuracy result, and the size gain it showed was largely an abundance-cutoff artifact (see Findings). |
| 4 | Label file has **several columns**; the first is primary. A blank or `NA` cell means unlabelled. The conflict check needs 3 or more members in a group. Each column gets its own colors, ARI and conflict check. | Both datasets needed more than one label (species/lineage, plus source). Labels enter after the distances, so extra columns are cheap. |
| 5 | **R1 and R2 are merged** per sample with the same grouping rules as `gmg run`; `--r1-only` is a fast mode. | Consistent with `run`; extra coverage reduces false low-coverage flags. |
| 6 | Input cap is **fixed bases per sample** (R1 and R2 combined), default about 1 Gb, documented as a default to adjust. The audit reports achieved coverage and warns when its spread across the cohort is large. | Equal reads did not give equal coverage (2x to 59x with a 6M-read cap). Target-coverage stopping is deferred. |
| 7 | Design and test for about **1,500 samples** with exact all-pairs; refuse above a configurable `--max-samples`, suggesting scaled 1000, a subsample or a split. | About 5e10 comparisons at n = 559 took minutes; n = 1,500 stays tractable; an inverted index is not needed yet. |
| 8 | **Hash-presence PCA is the default ordination**; PCoA of 1-ANI stays as `--ordination pcoa`, written to a file but not in the headline figure. | Hash PCA tracked SNP PCA better on all four datasets (table below). |
| 9 | **Clusters:** cut UPGMA at the largest merge-height gap only when gap1/gap2 is at least 2. Otherwise report "no single clear cut" and list the top 3 candidate k with their gaps. k = 1 is allowed. `--k` and `--cut` override. | The rule recovered every SNP-supported split and refused to invent more (table below). The threshold of 2 rests on four datasets; after the cutoff correction Bd's ratio is 1.96 (just under 2) and picks the k=5 that matches the coarse groups, so the threshold needs revisiting. |
| 10 | **Tree:** NJ on mash distance with unsupported branches collapsed to polytomies; unrooted by default, `--root ID` available. The collapse floor is the split-half noise floor (see below). | Within-lineage trees were combs; a fully resolved tree would imply resolution the data lacks. Midpoint rooting failed whenever an outlier was present. |
| 11 | **Resolution limit** comes from a split-half noise floor measured per run, used by both the NJ collapse and the resolution advisor. Experiment run 2026-11-04: the floor predicts which Bd pairs are indistinguishable when it is a high quantile (90th percentile), but not when it is the median, and the second check was inconclusive. See "Noise floor"; **needs revisiting before it is built**. Fallback is a fixed floor of about 10 differing hashes plus the gap ratio. | Bd nearest-neighbor ANI is about 0.9995 and lineages 2/3/4 cannot be told apart even though about 1% of hashes differ, so a hash-count floor is the wrong test. |
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
   Lecanicillium, NM_7898 Histoplasma). **Open: contamination.** With the proper cutoff a size test does not see contamination
   at about 1-20%, because contaminant k-mers fall below the cutoff (the 9 strains reported as mixtures by `gmg run` are all at
   normal size, 25.4-28.4 Mb). Detecting it needs another statistic (the abundance spectrum, or the index's existing
   `mixture_suspected` / `second_taxon`, which needs `--index`). This sits awkwardly with decision 2 (the index is never
   required) and needs a decision.
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

## Thresholds still to be tuned (not decided)

The MAD multiplier and floor for QC flags; the gap-ratio threshold of 2; the 1 Gb cap; the collapse floor and its rescaling;
the 30% flagged warning. Tune on the Cocci and Bd data; the fixture cannot set them.

## Build order (small PRs)

0. Two experiments, no product code. **Both run 2026-11-04.** (a) gmg-core's `auto_min_abund` on all 559 Cocci and 334 Bd
   scaled-200 sketches: it works (cutoffs 3-16 Cocci, mostly 4 Bd) and nothing needs porting; the Python prototype's rule
   returned 2 for almost every sample and was wrong. The 27 samples a first Python rule mishandled are low-coverage (2-9x), not
   high-coverage as first written; the Rust rule returns 2 for them too. (b) the split-half noise floor: see "Noise floor".
1. Prep PR: extract input discovery, grouping, sample-sheet parsing and resume logic into a small library crate used by
   `gmg run` and the new crate. Behavior-preserving; covered by the existing golden tests.
2. New crate core, with the simulated fixture written alongside so every later step has a test: QC, distances, UPGMA, NJ with
   collapse, Newick, clusters by the gap rule.
3. Hash PCA (and PCoA) with an eigensolver.
4. Labels: multi-column parser, ARI, conflict checks.
5. CLI and the output schema doc.
6. Report: scatter with hover, tree with colored tips, QC plot, tables, scalable heatmap.
7. Resolution advisor.
8. Validation on the HPCC (Cocci and Bd documented as regression checks) and a tutorial.

## Findings the design rests on

**Correction (2026-11-04).** The first explorations used a Python abundance cutoff that returned 2 for almost every sample.
Rerunning with gmg-core's `auto_min_abund` (no decontamination) changed several results, listed here. Numbers below are from
the reruns unless marked.

- Coccidioides: species gap in ANI (within about 0.996, between about 0.980); both species and the 19-sample immitis cluster
  are exact SNP clades (that test used the earlier cluster assignments); posadasii k=3 supported for one split only (12 of 262
  samples discordant); two `samples.csv` label problems (not re-checked with the corrected cutoff: the rerun applied the
  earlier relabelling).
- **Retracted:** "31 contaminated Coccidioides samples" and the decontamination gain (32.6 to 27.4 Mb). With the proper cutoff
  the cohort is tight (25.9 Mb median, 5-95% 25.5-26.7 at 8x or more); 25 of the 33 flagged samples are no longer outliers and
  the 8 that remain are mostly 2-5x samples. The 9 strains reported as mixtures are at normal size (25.4-28.4 Mb), so size
  cannot detect contamination at these levels.
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
  rule. Tested 2026-11-04: not borne out. The Rust rule works on all 559 Cocci and 334 Bd sketches; the Python rule is the
  one that is wrong (it returns 2 for almost every sample). The reviewer was right that the two rules differ.
- The suggested reuse of `classify.rs:12` (`READS_SPECIES_MIN_ANI`) does not apply: it lowers a species cutoff for
  classification and does not correct an ANI estimate.
- That min-size normalisation inflates the between-species ANI of 0.978 is a hypothesis to test, not a finding; the two
  species have similar genome sizes.
