# Headline precision re-measured on release 2026.10.5

Date 2026-10-31. Reference: release `index-2026.10.5` (`bfd-2026.10.5-s10000.gmgi`, 23,677 genomes = 23,669 BFD fungi minus the 72 ambiguous-label genomes of
`reports/label_audit/reference_exclusions.txt`, plus 80 non-fungal decoys; scaled = 10000, k = 21). Binary `gmg-ca8e3b9` (main at #119: genus vote consistency 0.9,
family ANI 0.7715, `mixture_suspected` and `high_error_reads` flags). Default (`known`) thresholds. Nothing was tuned; the scripts and the inputs are those of the
earlier reports. "Earlier" is release 2026.10.3 (truth set, with the binary of #89/#93 and the thresholds before #58) or 2026.10.x (BFD index with 43 exclusions,
2026.10.3 and earlier) as the cell says; the exclusions went from 43 (bench) or 57 (truth set, post-freeze) to 72, and the thresholds from genus vote consistency 0.5 to 0.9.

What was re-run (HPCC, `scripts/hpcc/40_*.sbatch`, plus the existing `23_truth_run.sbatch` with `TAG=i5` and `25_release_stability.sbatch` with `VERS=2026.10.5`):

| set | how | files here |
|---|---|---|
| RefSeq truth set, leave-one / species / genus-out (652 genomes, one per mode) and 21 outside the index | `23_truth_run.sbatch` (`gmg run --exclude-level self/species/genus`), scored with `validate_run.py` against the RefSeq lineage (`40_truth_score_i5.sbatch`); GCF_011074885.1 is one of the 72 exclusions, so it is not in the index and is dropped as before | `truth_set_tables.md`, `*.known.i5.raw.refseq.tsv` |
| Leave-out bench on the BFD genomes (`gmg bench`, loo / los / log, seed 1, all genomes) | `40_bench_i5.sbatch` (reference list = the 23,597 fungal genomes of the release, no decoys; the 72 excluded ids are dropped as queries), `novelty_curve.py`, `heldout_calibration.py` (20 splits, target 0.95) | `novelty_curve.*`, `heldout_calibration.*` |
| 1,257 post-freeze assemblies | `25_release_stability.sbatch`, `validate_run.py` | `post_freeze_balanced.*` |
| Real reads: 558 *Coccidioides* sets (R1, first 100k / 500k / 2M reads), 345 RNA-seq sets | `40_reads_i5.sbatch` (the commands of `11_cocci_reads` / `12_rna_reads`), `validate_run.py` | `cocci_reads_n*.md`, `rna_reads_n*.md` |

## Old against new, by rank

Precision of the calls made (called fraction in brackets). Moves of more than 1 point are marked **bold**.

**Truth set, leave-one-out (652)**

| | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|
| earlier (2026.10.3) | 0.978 (0.98) | 0.971 (0.98) | 0.969 (0.97) | 0.986 (0.92) | 0.994 (0.77) | 0.949 (0.39) |
| 2026.10.5 | 0.978 (0.98) | 0.971 (0.98) | 0.968 (0.97) | 0.988 (0.92) | 0.994 (0.77) | 0.949 (0.39) |

**Truth set, leave-species-out (652)**

| | phylum | class | order | family | genus |
|---|---|---|---|---|---|
| earlier | 0.962 (0.96) | 0.948 (0.97) | 0.944 (0.95) | 0.960 (0.87) | 0.980 (0.54) |
| 2026.10.5 | 0.962 (0.96) | 0.948 (0.97) | 0.944 (0.95) | 0.964 (0.86) | 0.983 (0.53) |

**Truth set, leave-genus-out (652)**

| | phylum | class | order | family | genus called |
|---|---|---|---|---|---|
| earlier | 0.860 (0.93) | 0.768 (0.94) | 0.707 (0.83) | 0.652 (0.57) | 47 (7%) |
| 2026.10.5 | 0.861 (0.93) | **0.783 (0.94)** | **0.726 (0.84)** | **0.687 (0.57)** | 47 (7%) |

The move at class (+1.5), order (+1.9) and family (+3.5) is attributed to the 15 further exclusions of 2026.10.5 (the class, order and family cutoffs and the decoy re-call are the same in both runs; no separate run isolates the exclusions). Wrong leave-genus-out calls: 523 against 554 disagreements, of which
`reference_label_error` 91 (17%) against 118 (21%). Leave-one-out and leave-species-out disagreement counts are unchanged (75; 158 against 161).
Outside the index (21): 21/21 correct at every rank with a truth, species called and correct for 18, as before.

**BFD leave-out bench (false calls for absent taxa; `novelty_curve.md`)**

| | earlier (43 exclusions) | 2026.10.5 (72 exclusions) |
|---|---|---|
| genus, ANI 0.90: leave-one-out coverage / precision | 0.906 / 0.985 | 0.906 / 0.986 |
| genus named when the genus is absent (leave-genus-out) | 0.183 | **0.105** |
| species, ANI 0.9938: leave-one-out coverage / precision | 0.722 / 0.982 | 0.724 / 0.982 |
| species named when the species is absent | 0.117 | 0.110 |

`novelty_curve.py` fixes the genus vote consistency at 0.5 (the value of its time); with the deployed 0.9 the numbers are coverage / precision 0.898 / 0.989 and
absent-genus false-call rate 0.102 (earlier bench 0.897 / 0.988 and 0.180). To check that the drop is the reference and not the script, the earlier bench was
re-scored with the 72 ids dropped as queries: 0.182, i.e. unchanged. The drop is therefore an effect of removing the 72 genomes as references (the genomes whose
labels conflicted with their relatives were being named for absent genera), measured on the same BFD labels. Caveat: the bench queries are BFD labels, which were not re-audited.

**Held-out calibration (`heldout_calibration.md`, `known` thresholds, 20 genus splits; precision / coverage on the held-out half)**

| rank | earlier | 2026.10.5 |
|---|---|---|
| phylum | 0.962 / 0.996 | 0.962 / 0.996 |
| class | 0.952 / 0.992 | 0.952 / 0.992 |
| order | 0.947 / 0.982 | 0.948 / 0.985 |
| family | 0.945 / 0.885 | 0.944 / **0.961** |
| genus | 0.946 / 0.978 | 0.946 / 0.980 |
| species | 0.948 / 0.548 | 0.949 / 0.539 |

(Clade-out precision, fitted on Ascomycota and tested on the other phyla: phylum 0.921, class 0.927, order 0.917, family 0.920 (was 0.945), genus 0.874, species 0.893.)
The family fit is the only held-out row that moved (coverage +7.6 points at the same precision, and clade-out family precision 0.945 to 0.920); the fit is a
calibration of the thresholds on the bench, not the shipped cutoffs.

**Independent assemblies (1,257 post-freeze, `post_freeze_balanced.md`)**

| | phylum | class | order | family | genus | species (52) |
|---|---|---|---|---|---|---|
| earlier | 0.990 | 0.991 | 0.986 | 0.983 | 0.988 | 1.000 (0.75) |
| 2026.10.5 | 0.990 | 0.991 | 0.988 | 0.984 | 0.988 | 1.000 (0.75) |

Called: phylum 0.946, genus 0.839 (both unchanged); `tentative_genus` inputs 78 to 74. Genus absent from the reference (22 genomes): genus named for 3, as before.

**Real reads (`cocci_reads_n*.md`, `rna_reads_n*.md`)**

| | 100k | 500k | 2M |
|---|---|---|---|
| *Coccidioides* genus called / precision, earlier | 0.996 / 0.995 | 1.000 / 0.993 | 0.998 / 0.995 |
| same, 2026.10.5 | 0.996 / 0.995 | 1.000 / 0.993 | 0.998 / 0.995 |
| *Coccidioides* species called / precision, earlier | 0.950 / 0.994 | 0.957 / 0.996 | 0.962 / 0.996 |
| same, 2026.10.5 | 0.950 / 0.994 | 0.957 / 0.996 | 0.962 / 0.996 |
| RNA-seq genus called / precision, earlier | 0.930 / 0.997 | 0.933 / 0.997 | 0.939 / 0.997 |
| same, 2026.10.5 | 0.924 / 0.997 | 0.924 / 1.000 | 0.930 / 1.000 |
| RNA-seq species called / precision, earlier | 0.678 / 0.966 | 0.684 / 0.966 | 0.629 / 0.968 |
| same, 2026.10.5 | 0.684 / 0.958 | 0.687 / 0.958 | 0.632 / 0.963 |

(Earlier values are the `*_reads099.md` reports. No cell moved by 1 point or more; RNA-seq genus called fell 0.6-0.9 points and species precision 0.5-0.8 points.)

## Not re-measured

The other limitations rows are unchanged by a reference swap of 72 genomes or were measured on inputs whose scoring needs data that is not part of this change
(hidden-species/genus RNA-seq, ZyGoLife and SRA read sets, mixtures, contamination and stress simulations, tool comparison, long reads, the cluster-based
mixture detector, IndigoAg and the other external sets). They keep their release label. The disagreement log `data/validation/disagreements.tsv` and the T2-T6
tables of `truth_set.md` (novel thresholds, decoy history, category tables) were not regenerated; the category counts on 2026.10.5 are in the text above.
