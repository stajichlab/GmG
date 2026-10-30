# Trusted truth set: NCBI RefSeq fungal genomes (#39)

Index `bfd-2026.10.3-s10000.gmgi` (23,612 BFD genomes without the 57 ambiguous-label ones, plus 81 non-fungal decoys, taxonomy pinned to NCBI
taxdump 2026-10-29), binary `gmg-a9acf5c` (decoy re-call #89, flag #93), default (`known`) thresholds unless noted, `gmg run --exclude-level` (#85).
Per-genome output stays on the HPCC (`GmG_runs/truthset/`); this page and `data/validation/` hold the aggregates, the truth list and the
disagreement log. The first run of this set (index 2026.10.2) is kept below as "before" where it explains a change.

## Summary

* **The truth set.** 673 RefSeq fungal assemblies: 267 from type material, 402 reference genomes, 4 other RefSeq; 68 chromosome, 68 complete,
  324 scaffold, 213 contig level. Truth = the NCBI lineage of the accession's taxid. 652 are in the index (the same assembly under its GenBank or RefSeq accession;
  258 are type material); 21 are not (9 type material) and are a plain, fully independent set. (One genome, GCF_011074885.1, was dropped from the set
  with the 57 exclusions.)
* **Leave-one-out (a new strain of a known species), 652 genomes, 258 of them type material.** Precision of the calls made, phylum to genus:
  0.978, 0.971, 0.969, 0.986, 0.994; species 0.949 (253 calls). For the 403 whose species has another genome in the reference the called
  fraction is 0.995 / 0.995 / 0.995 / 0.987 / 0.970 / 0.60 (phylum ... species) at precision 0.998 / 0.995 / 0.995 / 1.000 / 1.000 / 0.992.
  The wrong phylum to family calls are all `low_confidence` (14, 18, 19, 8 of them); of the 16 confident disagreements none is a plain tool error
  (13 species, 3 genus: sister species above the species cutoff, provisional "sp." labels, genus synonyms, one label conflict).
* **Leave-species-out (a new species of a known genus).** Precision phylum to genus 0.962 / 0.948 / 0.944 / 0.960 / 0.980. 42 inputs (6.4%) got a
  species call, necessarily a wrong one: 25 are a congener at ANI >= 0.99, 17 coincide with a label conflict or a provisional label; none was a
  plain overcall at lower ANI.
* **Leave-genus-out (a new genus of a known family).** Precision phylum 0.860, class 0.768, order 0.707, family 0.652 on the calls made (called
  0.93 / 0.94 / 0.83 / 0.57). 47 inputs (7.2%) got a genus call: 15 coincide with a label conflict, 20 are accepted synonym genera (*Aspergillus/Neosartorya*,
  *Botrytis/Botryotinia*, ...), 12 are plain overcalls at ANI 0.90-0.95.
* **Outside the index (21 assemblies).** Every call from phylum to genus is correct (21 of 21 at each rank with a truth); species called and correct for 18;
  the other 3 withheld (`low_agreement`: the near-best hits carry different species labels at ANI 0.9998).
* **Against the first run (index 2026.10.2).** Leaving the 57 genomes out and re-making distant decoy calls without the decoys together raise raw
  leave-genus-out class precision from 0.561 to 0.768, order from 0.557 to 0.707, family from 0.454 to 0.652, and phylum from 0.679 to 0.860 (no input is called
  *Chordata* any more). See the two sections on the 57 genomes and on the decoys for what each change did.

## Method

1. `scripts/hpcc/20_truth_summary.sbatch`: `datasets summary genome taxon fungi --assembly-source RefSeq` (674 records; `type_material`,
   `refseq_category`, level, N50). `21_truth_prep.sbatch`: truth lineage from the taxid with `normalize_taxonomy.py` and the taxdump of the
   index; link to the index genome of the same assembly by accession number (`scripts/make_truth_set.py`, list in
   `data/validation/truth_set.list.tsv`); sample sheets with the `ref` column (`make_truth_sheets.py`). `22_truth_fetch.sbatch`: download the genomes.
2. `23_truth_run.sbatch`: `gmg run --sample-sheet ... --exclude-level self|species|genus` (leave-one/species/genus-out; the genome, its species or
   its genus is hidden from the reference for that query) and a plain run for the 21 outside the index. Species and genus come from the
   lineage of the hidden genome in the index. None of the 652 genomes is one of the 57 left out of the reference. One input is in the index twice
   (GCF_001792695.1) and 14 leave-one-out inputs have an identical copy under another accession (`identical_to_reference`; scoring without them
   changes no number by more than 0.003).
3. `scripts/validate_run.py` for the tables (precision = correct among calls made; called = share of inputs with a call; Wilson 95% intervals
   below); `scripts/review_disagreements.py` for the disagreement log; `scripts/truth_set_tables.py` for the tables here.
4. The earlier "flag-aware" scoring (a phylum call that only names a distant decoy counted as no call, `scripts/blank_decoy_calls.py`) is no longer needed:
   on 2026.10.3 it blanks 0 of the 652 + 652 + 652 + 21 inputs, so raw and flag-aware tables are identical. The script stays for scoring older indexes.

## Results: all inputs

Cells are *called fraction / precision [95% interval] (number of calls)*.

**T1. All inputs, called fraction / precision [95% Wilson] (calls made), index 2026.10.3**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out | 652 | 0.98 / 0.978 [0.963, 0.987] (636) | 0.98 / 0.971 [0.955, 0.982] (626) | 0.97 / 0.969 [0.952, 0.980] (618) | 0.92 / 0.986 [0.973, 0.993] (586) | 0.77 / 0.994 [0.983, 0.998] (502) | 0.39 / 0.949 [0.914, 0.970] (253) |
| leave-species-out | 652 | 0.96 / 0.962 [0.943, 0.974] (624) | 0.97 / 0.948 [0.928, 0.963] (619) | 0.95 / 0.944 [0.922, 0.959] (604) | 0.87 / 0.960 [0.940, 0.973] (550) | 0.54 / 0.980 [0.959, 0.990] (349) | 0.06 / 0.000 [0.000, 0.084] (42) |
| leave-genus-out | 652 | 0.93 / 0.860 [0.830, 0.885] (605) | 0.94 / 0.768 [0.732, 0.800] (598) | 0.83 / 0.707 [0.666, 0.744] (525) | 0.57 / 0.652 [0.601, 0.699] (362) | 0.07 / 0.000 [0.000, 0.076] (47) | 0.00 / 0.000 [0.000, 0.562] (3) |
| outside the index (21) | 21 | 1.00 / 1.000 [0.845, 1.000] (21) | 1.00 / 1.000 [0.839, 1.000] (20) | 1.00 / 1.000 [0.839, 1.000] (20) | 1.00 / 1.000 [0.845, 1.000] (21) | 1.00 / 1.000 [0.845, 1.000] (21) | 0.86 / 1.000 [0.824, 1.000] (18) |


Subsets of leave-one-out:

**leave-one-out, species still in the reference (answerable)**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, species still in the reference (answerable) | 403 | 0.99 / 0.998 [0.986, 1.000] (401) | 0.99 / 0.995 [0.982, 0.999] (393) | 0.99 / 0.995 [0.982, 0.999] (392) | 0.99 / 1.000 [0.990, 1.000] (390) | 0.97 / 1.000 [0.990, 1.000] (391) | 0.60 / 0.992 [0.970, 0.998] (242) |

**leave-one-out, species a singleton in the reference (not answerable)**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, species a singleton in the reference (not answerable) | 249 | 0.94 / 0.945 [0.908, 0.967] (235) | 0.97 / 0.931 [0.891, 0.957] (233) | 0.94 / 0.925 [0.883, 0.953] (226) | 0.82 / 0.959 [0.922, 0.979] (196) | 0.45 / 0.973 [0.924, 0.991] (111) | 0.04 / 0.000 [0.000, 0.259] (11) |

**leave-one-out, type material only**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, type material only | 258 | 0.98 / 0.957 [0.924, 0.976] (253) | 0.99 / 0.940 [0.904, 0.963] (251) | 0.97 / 0.943 [0.907, 0.966] (247) | 0.88 / 0.973 [0.942, 0.987] (220) | 0.64 / 0.994 [0.966, 0.999] (165) | 0.32 / 0.928 [0.851, 0.966] (83) |

**leave-one-out, reference genomes (not type)**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, reference genomes (not type) | 390 | 0.97 / 0.992 [0.977, 0.997] (379) | 0.98 / 0.992 [0.976, 0.997] (371) | 0.97 / 0.986 [0.969, 0.994] (367) | 0.95 / 0.995 [0.980, 0.998] (362) | 0.85 / 0.994 [0.978, 0.998] (333) | 0.43 / 0.958 [0.917, 0.980] (168) |


Leave-species-out split by whether the genus is still in the reference:

**leave-species-out, genus still in the reference**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-species-out, genus still in the reference | 576 | 0.96 / 0.975 [0.958, 0.985] (554) | 0.98 / 0.967 [0.949, 0.979] (551) | 0.97 / 0.965 [0.946, 0.977] (543) | 0.91 / 0.979 [0.962, 0.988] (513) | 0.60 / 0.988 [0.971, 0.995] (346) | 0.07 / 0.000 [0.000, 0.084] (42) |

**leave-species-out, genus is a singleton (also absent)**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-species-out, genus is a singleton (also absent) | 76 | 0.92 / 0.857 [0.757, 0.921] (70) | 0.93 / 0.794 [0.684, 0.873] (68) | 0.84 / 0.754 [0.633, 0.845] (61) | 0.52 / 0.703 [0.542, 0.825] (37) | 0.04 / 0.000 [0.000, 0.562] (3) | 0.00 / n/a (0) |


`--thresholds novel` (queries that may be a new genus): genus, family and species are never called, as designed; what remains is more precise
at class and order and much less often called.

**T2. `--thresholds novel`, index 2026.10.3**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-species-out | 652 | 0.84 / 0.989 [0.976, 0.995] (548) | 0.53 / 0.994 [0.979, 0.998] (339) | 0.40 / 0.996 [0.978, 0.999] (255) | 0.00 / n/a (0) | 0.00 / n/a (0) | 0.00 / n/a (0) |
| leave-genus-out | 652 | 0.52 / 0.956 [0.928, 0.973] (339) | 0.10 / 0.778 [0.661, 0.863] (63) | 0.04 / 0.643 [0.458, 0.793] (28) | 0.00 / n/a (0) | 0.00 / n/a (0) | 0.00 / n/a (0) |
| outside the index | 21 | 1.00 / 1.000 [0.845, 1.000] (21) | 1.00 / 1.000 [0.839, 1.000] (20) | 1.00 / 1.000 [0.839, 1.000] (20) | 0.00 / n/a (0) | 0.00 / n/a (0) | 0.00 / n/a (0) |


Why species is not called for 161 of the 403 answerable leave-one-out inputs: 111 have a best ANI below the species cutoff (0.9938: a strain
of a known species but a genetically distant one), 39 have near-best hits that carry different species labels (`low_agreement`), 11 stop at a higher
rank for other reasons. A species call is right when made (precision 0.992).

The wrong phylum to order calls in leave-one-out (14 / 18 / 19) and, much more, in leave-genus-out are the vote of the nearest genomes of a
lineage the reference barely has; every one is `low_confidence`. Use `status` and the `tentative_genus` column, not the call alone, when the
reference may lack the lineage.

## Trusted truth against the BFD benchmark

The full-BFD leave-one-out benchmark on BFD labels (`reports/overnight_2026-10-27.md`, answerable queries, top-vote accuracy) reads phylum
0.985, class 0.979, order 0.974, family 0.970, genus 0.958, species 0.928. The trusted leave-one-out set, answerable queries, reads (recall,
calls and non-calls counted as in the tables) 0.993, 0.990, 0.990, 0.987, 0.970 and, for species, 0.596 with a withholding rule that keeps species
precision at 0.992. The two measure slightly different things (vote accuracy against calls with thresholds), but the trusted set is 1 to 2 points
better from phylum to family and equal at genus: **the BFD labels of the *query* genomes are not what limits the benchmark, and the label noise
among the references (now partly removed) cost it about that much**. For all 652 genomes in both sets the lineage of the BFD record is identical to the
NCBI lineage of the RefSeq accession at all six ranks (the index lineage is derived from the same taxid), so the gap on the query side is zero by construction.

The label errors that matter are on the reference side: a genome of another lineage carrying a wrong or uninformative name. The disagreement
log separates them: with the 57 left out, **13 of the 75 leave-one-out disagreements (17%), 43 of 161 leave-species-out ones (27%) and 118 of
554 leave-genus-out ones (21%) coincide with a label conflict** (ANI >= 0.95 to a genome of a different lineage at that rank, an audit-listed
neighbour, or a provisional "sp." label); before the exclusions the leave-genus-out share was 247 of 410 (60%). A further 2, 5 and 21 are accepted genus
synonyms or recent splits, and 7, 25 and 0 are sister species above the species cutoff. What is left, the plain tool errors, is 53, 88 and
415 rows, almost all `low_confidence` calls of a neighbouring taxon (see below).

## Disagreements

Every wrong call at a rank with a known truth is in `data/validation/disagreements.tsv` (838 rows over the three modes and the novel-threshold
runs; the column `prev_category` is the category the same row had in the first run, 2026.10.2). Categories: `tool_error`, `reference_label_error`,
`taxonomy_synonym`, `species_complex`, `low_quality_input`, `truth_label_error`, `contamination` (the last three were never assigned: no call fell below
N50 10 kb, and neither a truth-side error nor contamination can be established from sequence similarity and the lineage alone). How the rest were
decided: rule-based first (NCBI synonym of the two names; genus pairs I listed by hand in `review_disagreements.py`; ANI >= 0.95 across a different
lineage; audit-listed top hit; same-genus species at ANI >= 0.99), then by hand. **The individual verdicts of the first run
(`data/validation/truth_set_reviews.tsv`, keyed by genome, rank and call) carry over: 23 of the 75 leave-one-out rows, 43 of 161 leave-species-out and
29 of 554 leave-genus-out have one, and the others follow the rules.** The 52 leave-one-out rows without an individual verdict (one genus pair and 51 new ones)
are `low_confidence` calls of the phylum to family at best-hit ANI 0.74-0.85 that the first run had blanked as distant-decoy calls. I read them all:
they are chance best hits on a few large or odd reference genomes (*Puccinia polysora* and *Phakopsora pachyrhizi*, about 1 Gb rust genomes;
*Gigaspora margarita*; *Neocallimastix* sp. G1) and not labelling problems, except that the *Neocallimastix* G1 genome is on the size-outlier list, which
the audit-list rule turned into six `reference_label_error` rows. The 415 leave-genus-out `tool_error` rows are mostly such calls of the
sister family or class (Aspergillaceae/Trichocomaceae, Hypoxylaceae/Xylariaceae, Podosporaceae/Lasiosphaeriaceae, Suillaceae/Auriculariaceae ...);
some are taxonomic differences between sources, not errors, and a systematist should look before they are counted against the tool.

**T5. Disagreements (a wrong call at a rank with a known truth) by category**
| dataset | disagreements | tool_error | reference_label_error | taxonomy_synonym | species_complex | low_quality_input | truth_label_error | contamination | status ok (confident) | reviewed by hand |
|---|---|---|---|---|---|---|---|---|---|---|
| loo.known | 75 | 53 | 13 | 2 | 7 | 0 | 0 | 0 | 16 | 23 |
| los.known | 161 | 88 | 43 | 5 | 25 | 0 | 0 | 0 | 53 | 43 |
| log.known | 554 | 415 | 118 | 21 | 0 | 0 | 0 | 0 | 99 | 29 |
| los.novel | 9 | 4 | 5 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| log.novel | 39 | 17 | 22 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| indep.known | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |


**loo.known: by rank**
| rank | tool_error | reference_label_error | taxonomy_synonym | species_complex | low_quality_input | total | of which status ok |
|---|---|---|---|---|---|---|---|
| phylum | 12 | 2 | 0 | 0 | 0 | 14 | 0 |
| class | 16 | 2 | 0 | 0 | 0 | 18 | 0 |
| order | 17 | 2 | 0 | 0 | 0 | 19 | 0 |
| family | 8 | 0 | 0 | 0 | 0 | 8 | 0 |
| genus | 0 | 1 | 2 | 0 | 0 | 3 | 3 |
| species | 0 | 6 | 0 | 7 | 0 | 13 | 13 |

**los.known: by rank**
| rank | tool_error | reference_label_error | taxonomy_synonym | species_complex | low_quality_input | total | of which status ok |
|---|---|---|---|---|---|---|---|
| phylum | 18 | 6 | 0 | 0 | 0 | 24 | 0 |
| class | 25 | 7 | 0 | 0 | 0 | 32 | 1 |
| order | 27 | 7 | 0 | 0 | 0 | 34 | 1 |
| family | 17 | 4 | 1 | 0 | 0 | 22 | 2 |
| genus | 1 | 2 | 4 | 0 | 0 | 7 | 7 |
| species | 0 | 17 | 0 | 25 | 0 | 42 | 42 |

**log.known: by rank**
| rank | tool_error | reference_label_error | taxonomy_synonym | species_complex | low_quality_input | total | of which status ok |
|---|---|---|---|---|---|---|---|
| phylum | 65 | 20 | 0 | 0 | 0 | 85 | 0 |
| class | 110 | 29 | 0 | 0 | 0 | 139 | 13 |
| order | 127 | 27 | 0 | 0 | 0 | 154 | 16 |
| family | 101 | 24 | 1 | 0 | 0 | 126 | 20 |
| genus | 12 | 15 | 20 | 0 | 0 | 47 | 47 |
| species | 0 | 3 | 0 | 0 | 0 | 3 | 3 |


**What changed against the first log.**

**T6. Rows whose category changed against the 2026.10.2 log (known thresholds)**
Rows in both logs: 315 (311 with the same category). Rows only in the 2026.10.3 log: 475. Rows only in the 2026.10.2 log (no longer a disagreement): 210 (reference_label_error 200, tool_error 10).

| 2026.10.2 | 2026.10.3 | rows |
|---|---|---|
| (not a disagreement in 2026.10.2) | tool_error | 396 |
| (not a disagreement in 2026.10.2) | reference_label_error | 77 |
| tool_error | reference_label_error | 3 |
| (not a disagreement in 2026.10.2) | taxonomy_synonym | 2 |
| reference_label_error | tool_error | 1 |


The 210 rows that disappeared are almost all `reference_label_error` rows whose neighbour is one of the 57 genomes now left out; the 475 new rows
are calls the first run had hidden (the blanked decoy phylum calls and their lineage) or made differently: with the decoy re-call the phylum to
order vote comes from fungal genomes again, and it is wrong more often than the blank was (see the decoy section).

## Effect of leaving the 57 ambiguous-label genomes out

Release 2026.10.3 leaves out 57 genomes (`reports/label_audit/reference_exclusions.txt`: the 43 of the label audit and the 14 this truth set found).
Raw tables, same inputs, same thresholds. "2026.10.2" is the first run (decoy phylum calls counted); "with the 43 hidden" is the first run with
`--exclude-ids`; "2026.10.3" is the production index with the decoy re-call:

**T4. Before and after leaving the 57 ambiguous-label genomes out**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, 2026.10.2 (43 audit genomes kept) | 653 | 0.98 / 0.950 [0.930, 0.964] (642) | 0.99 / 0.952 [0.933, 0.966] (629) | 0.97 / 0.961 [0.943, 0.974] (618) | 0.92 / 0.971 [0.954, 0.982] (587) | 0.77 / 0.992 [0.980, 0.997] (500) | 0.39 / 0.948 [0.914, 0.970] (252) |
| leave-one-out, 2026.10.2 with the 43 hidden | 653 | 0.98 / 0.950 [0.930, 0.964] (642) | 0.99 / 0.954 [0.935, 0.968] (629) | 0.97 / 0.961 [0.943, 0.974] (619) | 0.93 / 0.971 [0.954, 0.982] (591) | 0.77 / 0.992 [0.980, 0.997] (503) | 0.39 / 0.949 [0.914, 0.970] (253) |
| leave-one-out, 2026.10.3 (57 left out) | 652 | 0.98 / 0.978 [0.963, 0.987] (636) | 0.98 / 0.971 [0.955, 0.982] (626) | 0.97 / 0.969 [0.952, 0.980] (618) | 0.92 / 0.986 [0.973, 0.993] (586) | 0.77 / 0.994 [0.983, 0.998] (502) | 0.39 / 0.949 [0.914, 0.970] (253) |
| leave-species-out, 2026.10.2 (43 audit genomes kept) | 653 | 0.97 / 0.905 [0.880, 0.926] (633) | 0.98 / 0.905 [0.880, 0.926] (622) | 0.94 / 0.923 [0.899, 0.942] (597) | 0.87 / 0.924 [0.899, 0.943] (553) | 0.53 / 0.962 [0.937, 0.978] (346) | 0.07 / 0.000 [0.000, 0.079] (45) |
| leave-species-out, 2026.10.2 with the 43 hidden | 653 | 0.97 / 0.904 [0.878, 0.924] (633) | 0.98 / 0.908 [0.883, 0.929] (622) | 0.94 / 0.923 [0.899, 0.942] (600) | 0.88 / 0.925 [0.900, 0.944] (559) | 0.54 / 0.966 [0.941, 0.980] (350) | 0.07 / 0.000 [0.000, 0.080] (44) |
| leave-species-out, 2026.10.3 (57 left out) | 652 | 0.96 / 0.962 [0.943, 0.974] (624) | 0.97 / 0.948 [0.928, 0.963] (619) | 0.95 / 0.944 [0.922, 0.959] (604) | 0.87 / 0.960 [0.940, 0.973] (550) | 0.54 / 0.980 [0.959, 0.990] (349) | 0.06 / 0.000 [0.000, 0.084] (42) |
| leave-genus-out, 2026.10.2 (43 audit genomes kept) | 653 | 0.95 / 0.679 [0.641, 0.715] (620) | 0.95 / 0.561 [0.521, 0.600] (606) | 0.85 / 0.557 [0.515, 0.599] (540) | 0.68 / 0.454 [0.408, 0.501] (434) | 0.12 / 0.000 [0.000, 0.047] (78) | 0.01 / 0.000 [0.000, 0.299] (9) |
| leave-genus-out, 2026.10.2 with the 43 hidden | 653 | 0.95 / 0.669 [0.631, 0.705] (620) | 0.95 / 0.635 [0.596, 0.672] (605) | 0.85 / 0.600 [0.559, 0.641] (543) | 0.68 / 0.528 [0.481, 0.574] (434) | 0.10 / 0.000 [0.000, 0.054] (67) | 0.01 / 0.000 [0.000, 0.354] (7) |
| leave-genus-out, 2026.10.3 (57 left out) | 652 | 0.93 / 0.860 [0.830, 0.885] (605) | 0.94 / 0.768 [0.732, 0.800] (598) | 0.83 / 0.707 [0.666, 0.744] (525) | 0.57 / 0.652 [0.601, 0.699] (362) | 0.07 / 0.000 [0.000, 0.076] (47) | 0.00 / 0.000 [0.000, 0.562] (3) |


Reading the table. For leave-genus-out the production numbers are class 0.768, order 0.707, family 0.652 (phylum 0.860) at called fractions 0.94 / 0.83 /
0.57. Two changes act together between "2026.10.2" and "2026.10.3": the 57 exclusions and the decoy re-call. Index 2026.10.1 has no decoys and still
holds all 57 (table in the next section): its leave-genus-out numbers 0.861 / 0.686 / 0.659 / 0.573 against 0.860 / 0.768 / 0.707 / 0.652 isolate the
exclusions, **+8 points of class, +5 of order, +8 of family**, with the same phylum precision. Leave-one-out and leave-species-out move by at most 0.02 (genus
leave-species-out 0.962 to 0.980); the number of leave-genus-out genus calls falls from 78 to 47 and the label-conflict rows from 247 to 118.
The "about 0.90" class precision quoted before (docs/limitations.md) was the flag-aware figure for the 43 hidden *at a called fraction of 0.67*; the
production number at the called fraction of 0.94 is 0.768.

## Effect of the decoy panel (#79)

Index 2026.10.1 (no decoys) against 2026.10.2 (81 decoys, raw and flag-aware) and 2026.10.3 (decoys, call re-made without them when the nearest
genomes are decoys at ANI < 0.90, #89); phylum to genus, same taxonomy:

**T3. Decoys: phylum to genus, same taxonomy**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, 2026.10.1 (no decoys) | 653 | 0.98 / 0.978 [0.964, 0.987] (638) | 0.99 / 0.968 [0.951, 0.979] (628) | 0.97 / 0.968 [0.951, 0.979] (619) | 0.92 / 0.985 [0.971, 0.992] (584) | 0.77 / 0.992 [0.980, 0.997] (500) | 0.39 / 0.948 [0.914, 0.970] (252) |
| leave-one-out, 2026.10.2 (decoys, raw) | 653 | 0.98 / 0.950 [0.930, 0.964] (642) | 0.99 / 0.952 [0.933, 0.966] (629) | 0.97 / 0.961 [0.943, 0.974] (618) | 0.92 / 0.971 [0.954, 0.982] (587) | 0.77 / 0.992 [0.980, 0.997] (500) | 0.39 / 0.948 [0.914, 0.970] (252) |
| leave-one-out, 2026.10.2 (decoys, flag-aware) | 653 | 0.94 / 0.998 [0.991, 1.000] (611) | 0.95 / 0.993 [0.983, 0.997] (603) | 0.94 / 0.995 [0.985, 0.998] (597) | 0.90 / 0.997 [0.987, 0.999] (572) | 0.77 / 0.992 [0.980, 0.997] (500) | 0.39 / 0.948 [0.914, 0.970] (252) |
| leave-one-out, 2026.10.3 (decoys, call re-made without them) | 652 | 0.98 / 0.978 [0.963, 0.987] (636) | 0.98 / 0.971 [0.955, 0.982] (626) | 0.97 / 0.969 [0.952, 0.980] (618) | 0.92 / 0.986 [0.973, 0.993] (586) | 0.77 / 0.994 [0.983, 0.998] (502) | 0.39 / 0.949 [0.914, 0.970] (253) |
| leave-species-out, 2026.10.1 (no decoys) | 653 | 0.96 / 0.962 [0.944, 0.974] (626) | 0.97 / 0.942 [0.921, 0.958] (621) | 0.95 / 0.939 [0.917, 0.955] (606) | 0.87 / 0.951 [0.930, 0.966] (551) | 0.53 / 0.962 [0.937, 0.978] (346) | 0.07 / 0.000 [0.000, 0.079] (45) |
| leave-species-out, 2026.10.2 (decoys, raw) | 653 | 0.97 / 0.905 [0.880, 0.926] (633) | 0.98 / 0.905 [0.880, 0.926] (622) | 0.94 / 0.923 [0.899, 0.942] (597) | 0.87 / 0.924 [0.899, 0.943] (553) | 0.53 / 0.962 [0.937, 0.978] (346) | 0.07 / 0.000 [0.000, 0.079] (45) |
| leave-species-out, 2026.10.2 (decoys, flag-aware) | 653 | 0.88 / 0.997 [0.987, 0.999] (575) | 0.90 / 0.983 [0.968, 0.990] (573) | 0.88 / 0.984 [0.970, 0.992] (560) | 0.82 / 0.977 [0.960, 0.987] (523) | 0.53 / 0.962 [0.937, 0.978] (346) | 0.07 / 0.000 [0.000, 0.079] (45) |
| leave-species-out, 2026.10.3 (decoys, call re-made without them) | 652 | 0.96 / 0.962 [0.943, 0.974] (624) | 0.97 / 0.948 [0.928, 0.963] (619) | 0.95 / 0.944 [0.922, 0.959] (604) | 0.87 / 0.960 [0.940, 0.973] (550) | 0.54 / 0.980 [0.959, 0.990] (349) | 0.06 / 0.000 [0.000, 0.084] (42) |
| leave-genus-out, 2026.10.1 (no decoys) | 653 | 0.94 / 0.861 [0.831, 0.886] (612) | 0.95 / 0.686 [0.648, 0.722] (602) | 0.86 / 0.659 [0.619, 0.698] (546) | 0.61 / 0.573 [0.523, 0.621] (386) | 0.12 / 0.000 [0.000, 0.047] (78) | 0.01 / 0.000 [0.000, 0.299] (9) |
| leave-genus-out, 2026.10.2 (decoys, raw) | 653 | 0.95 / 0.679 [0.641, 0.715] (620) | 0.95 / 0.561 [0.521, 0.600] (606) | 0.85 / 0.557 [0.515, 0.599] (540) | 0.68 / 0.454 [0.408, 0.501] (434) | 0.12 / 0.000 [0.000, 0.047] (78) | 0.01 / 0.000 [0.000, 0.299] (9) |
| leave-genus-out, 2026.10.2 (decoys, flag-aware) | 653 | 0.67 / 0.966 [0.944, 0.979] (436) | 0.68 / 0.785 [0.744, 0.821] (433) | 0.62 / 0.766 [0.722, 0.805] (393) | 0.50 / 0.616 [0.561, 0.667] (320) | 0.12 / 0.000 [0.000, 0.047] (78) | 0.01 / 0.000 [0.000, 0.299] (9) |
| leave-genus-out, 2026.10.3 (decoys, call re-made without them) | 652 | 0.93 / 0.860 [0.830, 0.885] (605) | 0.94 / 0.768 [0.732, 0.800] (598) | 0.83 / 0.707 [0.666, 0.744] (525) | 0.57 / 0.652 [0.601, 0.699] (362) | 0.07 / 0.000 [0.000, 0.076] (47) | 0.00 / 0.000 [0.000, 0.562] (3) |


On 2026.10.2 a fungus with no close relative had a decoy as nearest genome by chance (*Mus musculus*, *Thalia democratica*, *Lampetra planeri*; ANI
0.75-0.88, a few dozen shared hashes) and its phylum vote read "Chordata": 31 of 653 inputs in leave-one-out, 58 in leave-species-out and 184 (28%) in
leave-genus-out, all `low_confidence`. Raw phylum precision fell from 0.978 to 0.950, 0.962 to 0.905 and 0.861 to 0.679. The flag-aware scoring
(a distant decoy call counted as no call) restored precision (0.998 / 0.997 / 0.966) by giving up the call (called 0.94 / 0.88 / 0.67). **With the
re-call on 2026.10.3 the numbers are those of the index without decoys again: phylum 0.978 / 0.962 / 0.860 at called 0.98 / 0.96 / 0.93.** No input is called
a decoy phylum in any mode (0 of 652 in each), and none of the fungi is flagged `not_fungal`. The call made without the decoys is the same low-confidence guess
the no-decoy index makes: right for about half to 60% of those inputs, so the blank gave a higher precision at a lower called fraction and the re-call gives a higher called fraction at
the precision of the no-decoy index; the choice is a policy one, and both are consistent with `status = low_confidence`. The flag `distant_non_fungal_hit` is still
set on 32 (4.9%), 60 (9.2%) and 202 (31%) inputs: the nearest genome is a distant decoy, a phylum call was made, and it came from the fungal genomes.

**The mechanism is not specific to decoys.** The 85 wrong leave-genus-out phylum calls (24 leave-species-out, 14 leave-one-out) have best-hit ANI 0.74-0.86
(median 0.79; 82 of 85 below 0.85), and their best hits are concentrated on a few genomes: *Neocallimastix californiae* (22) and *N. lanati* (14),
*Gigaspora margarita* (11), *Puccinia polysora* (7), *Massospora platypediae* (6). A genome with a very large or repeat-rich assembly shares a few dozen
hashes with almost any divergent query by chance, whether it is a mouse or a rust fungus. By best-hit ANI the phylum call has this precision (calls made, n in brackets):

| best-hit ANI | leave-one-out | leave-species-out | leave-genus-out |
|---|---|---|---|
| < 0.80 | 0.852 (54) | 0.854 (82) | 0.817 (279) |
| 0.80-0.85 | 0.857 (35) | 0.889 (90) | 0.818 (170) |
| 0.85-0.90 | 0.977 (43) | 0.978 (92) | 0.964 (83) |
| 0.90-0.95 | 1.000 (57) | 1.000 (109) | 1.000 (46) |
| >= 0.95 | 1.000 (447) | 1.000 (251) | 1.000 (27) |

So a phylum call is trustworthy from best-hit ANI 0.85 up (precision 0.96-1.00) and a guess below it (0.82-0.89); `status` already says `low_confidence` there,
and the results table could say it in the call itself (suggestion for a later change: a "phylum, weak" marker, or no phylum below 0.85).

## Outside the index

The 21 RefSeq assemblies whose accession number is not in the index (model organisms such as *Pyricularia oryzae*, *Candida albicans*,
*Trichoderma reesei*: the species is in the index through another assembly). All 21 `ok`; correct at phylum, class, order, family and
genus wherever a truth exists; species called and correct for 18 of 21 (first run: 17); the 3 withheld are *Fusarium verticillioides*,
*Paracoccidioides brasiliensis* and *Trichoderma virens* at ANI 0.9998-0.9999 with `low_agreement` (near-best hits of different species labels:
label inconsistency in the reference, not a failure of the match). Nine are type material.

## Limits

* The truth is the NCBI lineage as of the pinned taxdump; a taxonomy change that NCBI has not made is counted against the tool
  (*Pichiomycetes*, genus splits). Where a name pair is a known synonym or split I categorised it `taxonomy_synonym` (list in
  `review_disagreements.py`, my reading, not an authority).
* The truth set is dominated by well-studied genera; only 21 of 673 are outside the index, so the independent estimate is small.
* "Leave-one-out" hides one accession; identical copies under other accessions remain (14 inputs, no material effect). It is not
  leave-strain-out.
* The leave-genus-out rule hides the genus of the genome as the index holds it; a lineage whose genus is polyphyletic in NCBI keeps relatives.
* Label conflicts are a similarity argument: two genomes at ANI >= 0.95 cannot both be right about a different genus, but which of the two is
  wrong, or whether the assembly is contaminated, was not checked against the assembly. Report, do not delete.
* The 52 leave-one-out rows and most of the 415 leave-genus-out `tool_error` rows have rule verdicts, not individual ones; the tables' counts are exact, the verdicts are not.
* "Before" numbers are the first run on index 2026.10.2 with the binary of that date; the 2026.10.1 numbers are from the first run too.

## Acceptance (#39)

| criterion | status |
|---|---|
| >= 250 type-material genomes scored in all three modes | 258 in leave-one-out, leave-species-out and leave-genus-out |
| every disagreement categorised | yes: 838 rows (known and novel runs), 4 categories used; 95 of the 790 known-threshold rows carry an individual verdict, the rest rule verdicts |
| share of BFD-benchmark errors that were label errors | query side 0 (labels identical); reference side 17% / 27% / 21% of the disagreements with the 57 left out (37% / 45% / 60% before) |
| accuracy on the trusted set next to the full-BFD benchmark | above |
