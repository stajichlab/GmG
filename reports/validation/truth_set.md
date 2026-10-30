# Trusted truth set: NCBI RefSeq fungal genomes (#39)

Index `bfd-2026.10.2-s10000.gmgi` (23,669 BFD genomes + 81 non-fungal decoys, taxonomy pinned to NCBI taxdump 2026-10-29), default
(`known`) thresholds unless noted, `gmg run --exclude-level` (#85). Per-genome output stays on the HPCC
(`GmG_runs/truthset/`); this page and `data/validation/` hold the aggregates, the truth list and the disagreement log.

## Summary

* **The truth set.** All 674 RefSeq fungal assemblies: 267 from type material (258 are in the index), 668 reference/representative genomes;
  68 chromosome, 69 complete, 324 scaffold, 213 contig level. Truth = the NCBI lineage of the accession's taxid. 653 are in the index
  (the same assembly under its GenBank or RefSeq accession); 21 are not (9 of them type material) and are a plain, fully independent set.
* **Leave-one-out (a new strain of a known species), 653 genomes, 258 of them type material.** Precision of the calls made, phylum to genus:
  0.998, 0.993, 0.995, 0.997, 0.992; species 0.948 (252 calls). For the 404 whose species has another genome in the reference the called
  fraction is 0.99 / 0.99 / 0.99 / 0.98 / 0.96 / 0.60 (phylum ... species) at precision 1.000 / 0.995 / 0.997 / 0.997 / 0.997 / 0.992.
  All 8 plain tool errors (phylum to family) were `low_confidence`; of the 18 confident disagreements none was a plain tool error
  (13 species, 4 genus, 1 family: sister species above the species cutoff, provisional "sp." labels, genus synonyms and label conflicts).
* **Leave-species-out (a new species of a known genus).** Precision phylum to genus 0.997 / 0.983 / 0.984 / 0.977 / 0.962. 45 inputs (6.9%) got a
  species call, necessarily a wrong one: 25 are a congener at ANI >= 0.99, 20 coincide with a label conflict or a provisional label; none was a
  plain overcall at lower ANI.
* **Leave-genus-out (a new genus of a known family).** Precision phylum 0.966, class 0.785, order 0.766, family 0.616 on the calls made (called
  0.67 / 0.68 / 0.62 / 0.50). 78 inputs (12%) got a genus call: 51 coincide with a label conflict (ANI >= 0.95 to a genome of another genus or an
  audit-listed neighbour), 18 are accepted synonym genera (*Aspergillus/Neosartorya*, *Botrytis/Botryotinia*, ...), 9 are plain overcalls at ANI 0.90-0.94.
* **Outside the index (21 assemblies).** Every call from phylum to genus is correct (21 of 21 at each rank with a truth); species called and correct
  for 17; the other 4 withheld (`low_agreement`: the near-best hits carry different species labels at ANI 0.9997).
* **Two findings that are not about the accuracy of the classifier core.** (1) The decoy panel turns distant fungi into "Chordata" calls (4.7% of
  leave-one-out inputs, 28% of leave-genus-out inputs): see *Effect of the decoy panel*. (2) The 43 suspect genomes kept in the production index cost
  11 points of class precision in leave-genus-out (0.785 to 0.897 when hidden): see *Effect of the 43 suspect reference genomes*.

## Method

1. `scripts/hpcc/20_truth_summary.sbatch`: `datasets summary genome taxon fungi --assembly-source RefSeq` (674 records; `type_material`,
   `refseq_category`, level, N50). `21_truth_prep.sbatch`: truth lineage from the taxid with `normalize_taxonomy.py` and the taxdump of the
   index; link to the index genome of the same assembly by accession number (`scripts/make_truth_set.py`, list in
   `data/validation/truth_set.list.tsv`); sample sheets with the `ref` column (`make_truth_sheets.py`). `22_truth_fetch.sbatch`: download the 674 genomes.
2. `23_truth_run.sbatch`: `gmg run --sample-sheet ... --exclude-level self|species|genus` (leave-one/species/genus-out; the genome, its species or
   its genus is hidden from the reference for that query) and a plain run for the 21 outside the index. Species and genus come from the
   lineage of the hidden genome in the index. One input is in the index twice (GCF_001792695.1) and 14 leave-one-out inputs have an identical copy
   under another accession (`identical_to_reference`; scoring without them changes no number by more than 0.003).
3. `scripts/validate_run.py` for the tables (precision = correct among calls made; called = share of inputs with a call; Wilson 95% intervals
   below); `scripts/review_disagreements.py` for the disagreement log; `scripts/truth_set_tables.py` for the tables here.
4. **Flag-aware scoring.** A phylum call that only says "nearest to a non-fungal decoy at ANI 0.75-0.88" is, by the flag's own documentation
   (`distant_non_fungal_hit`, docs/qc.md), not a verdict. The "flag-aware" tables treat those calls as no call (`scripts/blank_decoy_calls.py`);
   the "raw" tables count them. Both are shown.

## Results: all inputs

Cells are *called fraction / precision [95% interval] (number of calls)*.

**T1. All inputs, called fraction / precision [95% Wilson] (calls made), flag-aware**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out | 653 | 0.94 / 0.998 [0.991, 1.000] (611) | 0.95 / 0.993 [0.983, 0.997] (603) | 0.94 / 0.995 [0.985, 0.998] (597) | 0.90 / 0.997 [0.987, 0.999] (572) | 0.77 / 0.992 [0.980, 0.997] (500) | 0.39 / 0.948 [0.914, 0.970] (252) |
| leave-species-out | 653 | 0.88 / 0.997 [0.987, 0.999] (575) | 0.90 / 0.983 [0.968, 0.990] (573) | 0.88 / 0.984 [0.970, 0.992] (560) | 0.82 / 0.977 [0.960, 0.987] (523) | 0.53 / 0.962 [0.937, 0.978] (346) | 0.07 / 0.000 [0.000, 0.079] (45) |
| leave-genus-out | 653 | 0.67 / 0.966 [0.944, 0.979] (436) | 0.68 / 0.785 [0.744, 0.821] (433) | 0.62 / 0.766 [0.722, 0.805] (393) | 0.50 / 0.616 [0.561, 0.667] (320) | 0.12 / 0.000 [0.000, 0.047] (78) | 0.01 / 0.000 [0.000, 0.299] (9) |
| outside the index (21) | 21 | 1.00 / 1.000 [0.845, 1.000] (21) | 1.00 / 1.000 [0.839, 1.000] (20) | 1.00 / 1.000 [0.839, 1.000] (20) | 1.00 / 1.000 [0.845, 1.000] (21) | 1.00 / 1.000 [0.845, 1.000] (21) | 0.81 / 1.000 [0.816, 1.000] (17) |


Raw (the decoy phylum call counted as a call):

**T1b. Same, raw (a distant decoy phylum call counted as a call)**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out | 653 | 0.98 / 0.950 [0.930, 0.964] (642) | 0.99 / 0.952 [0.933, 0.966] (629) | 0.97 / 0.961 [0.943, 0.974] (618) | 0.92 / 0.971 [0.954, 0.982] (587) | 0.77 / 0.992 [0.980, 0.997] (500) | 0.39 / 0.948 [0.914, 0.970] (252) |
| leave-species-out | 653 | 0.97 / 0.905 [0.880, 0.926] (633) | 0.98 / 0.905 [0.880, 0.926] (622) | 0.94 / 0.923 [0.899, 0.942] (597) | 0.87 / 0.924 [0.899, 0.943] (553) | 0.53 / 0.962 [0.937, 0.978] (346) | 0.07 / 0.000 [0.000, 0.079] (45) |
| leave-genus-out | 653 | 0.95 / 0.679 [0.641, 0.715] (620) | 0.95 / 0.561 [0.521, 0.600] (606) | 0.85 / 0.557 [0.515, 0.599] (540) | 0.68 / 0.454 [0.408, 0.501] (434) | 0.12 / 0.000 [0.000, 0.047] (78) | 0.01 / 0.000 [0.000, 0.299] (9) |
| outside the index (21) | 21 | 1.00 / 1.000 [0.845, 1.000] (21) | 1.00 / 1.000 [0.839, 1.000] (20) | 1.00 / 1.000 [0.839, 1.000] (20) | 1.00 / 1.000 [0.845, 1.000] (21) | 1.00 / 1.000 [0.845, 1.000] (21) | 0.81 / 1.000 [0.816, 1.000] (17) |


Subsets of leave-one-out:

**leave-one-out, species still in the reference (answerable)**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, species still in the reference (answerable) | 404 | 0.99 / 1.000 [0.990, 1.000] (400) | 0.99 / 0.995 [0.982, 0.999] (392) | 0.99 / 0.997 [0.986, 1.000] (390) | 0.98 / 0.997 [0.986, 1.000] (387) | 0.96 / 0.997 [0.986, 1.000] (389) | 0.60 / 0.992 [0.970, 0.998] (241) |

**leave-one-out, species a singleton in the reference (not answerable)**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, species a singleton in the reference (not answerable) | 249 | 0.85 / 0.995 [0.974, 0.999] (211) | 0.88 / 0.991 [0.966, 0.997] (211) | 0.86 / 0.990 [0.965, 0.997] (207) | 0.77 / 0.995 [0.970, 0.999] (185) | 0.45 / 0.973 [0.924, 0.991] (111) | 0.04 / 0.000 [0.000, 0.259] (11) |

**leave-one-out, type material only**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, type material only | 258 | 0.90 / 1.000 [0.984, 1.000] (233) | 0.91 / 0.991 [0.969, 0.998] (232) | 0.91 / 0.991 [0.969, 0.998] (231) | 0.84 / 1.000 [0.982, 1.000] (210) | 0.64 / 0.994 [0.966, 0.999] (165) | 0.32 / 0.928 [0.851, 0.966] (83) |

**leave-one-out, reference genomes (not type)**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, reference genomes (not type) | 391 | 0.96 / 0.997 [0.985, 1.000] (374) | 0.97 / 0.995 [0.980, 0.999] (367) | 0.96 / 0.997 [0.985, 1.000] (362) | 0.94 / 0.994 [0.980, 0.998] (358) | 0.85 / 0.991 [0.974, 0.997] (331) | 0.43 / 0.958 [0.916, 0.980] (167) |


Leave-species-out split by whether the genus is still in the reference:

**leave-species-out, genus still in the reference**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-species-out, genus still in the reference | 577 | 0.92 / 0.996 [0.986, 0.999] (532) | 0.94 / 0.981 [0.966, 0.990] (531) | 0.92 / 0.985 [0.970, 0.992] (520) | 0.88 / 0.986 [0.971, 0.993] (495) | 0.59 / 0.974 [0.951, 0.986] (342) | 0.08 / 0.000 [0.000, 0.079] (45) |

**leave-species-out, genus is a singleton (also absent)**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-species-out, genus is a singleton (also absent) | 76 | 0.57 / 1.000 [0.918, 1.000] (43) | 0.58 / 1.000 [0.916, 1.000] (42) | 0.55 / 0.975 [0.871, 0.996] (40) | 0.39 / 0.821 [0.644, 0.921] (28) | 0.05 / 0.000 [0.000, 0.490] (4) | 0.00 / n/a (0) |


`--thresholds novel` (queries that may be a new genus): genus, family and species are never called, as designed; what remains is more precise
at the high ranks and much less often called.

**T2. `--thresholds novel`, flag-aware**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-species-out | 653 | 0.83 / 1.000 [0.993, 1.000] (540) | 0.54 / 0.986 [0.966, 0.994] (344) | 0.41 / 0.985 [0.961, 0.994] (262) | 0.00 / n/a (0) | 0.00 / n/a (0) | 0.00 / n/a (0) |
| leave-genus-out | 653 | 0.51 / 0.994 [0.978, 0.998] (335) | 0.15 / 0.628 [0.527, 0.719] (94) | 0.09 / 0.525 [0.400, 0.647] (59) | 0.00 / n/a (0) | 0.00 / n/a (0) | 0.00 / n/a (0) |
| outside the index | 21 | 1.00 / 1.000 [0.845, 1.000] (21) | 1.00 / 1.000 [0.839, 1.000] (20) | 1.00 / 1.000 [0.839, 1.000] (20) | 0.00 / n/a (0) | 0.00 / n/a (0) | 0.00 / n/a (0) |


Why species is not called for 163 of the 404 answerable leave-one-out inputs: 110 have a best ANI below the species cutoff (0.9938: a strain
of a known species but a genetically distant one), 40 have near-best hits that carry different species labels (`low_agreement`), 13 stop at a higher
rank for other reasons. A species call is right when made (precision 0.992).

## Trusted truth against the BFD benchmark

The full-BFD leave-one-out benchmark on BFD labels (`reports/overnight_2026-10-27.md`, answerable queries, top-vote accuracy) reads phylum
0.985, class 0.979, order 0.974, family 0.970, genus 0.958, species 0.928. The trusted leave-one-out set, answerable queries, reads (recall,
calls and non-calls counted as in the table above) 0.990, 0.985, 0.985, 0.975, 0.960 and, for species, 0.592 with a withholding rule that keeps
species precision at 0.992. The two measure slightly different things (vote accuracy against calls with thresholds), but they are within a
point of each other from phylum to genus: **the BFD labels of the *query* genomes are not what limits the benchmark**. This is not an accident
of the comparison: for all 653 genomes in both sets the lineage of the BFD record is identical to the NCBI lineage of the RefSeq accession at all
six ranks (the index lineage is derived from the same taxid), so the gap on the query side is zero by construction.

The label errors that matter are on the reference side: a genome of another lineage carrying a wrong or uninformative name. The disagreement
log separates them (next section): **37% of the leave-one-out disagreements (10 of 27), 45% of the leave-species-out ones (41 of 91) and 60% of the
leave-genus-out ones (247 of 410) coincide with a label conflict** (ANI >= 0.95 to a genome of a different lineage at that rank, an audit-listed neighbour, or a
provisional "sp." label); a further 2, 5 and 19 are accepted genus synonyms or recent splits, and 7, 25 and 0 are sister species above the species
cutoff. What is left, the plain tool errors, is 8, 20 and 144.

## Disagreements

Every wrong call at a rank with a known truth is in `data/validation/disagreements.tsv` (602 rows over the three modes and the novel-threshold
runs). Categories: `tool_error`, `reference_label_error`, `taxonomy_synonym`, `species_complex`, `low_quality_input`, `truth_label_error`,
`contamination` (the last three were never assigned: no call fell below N50 10 kb, and neither a truth-side error nor contamination can be
established from sequence similarity and the lineage alone). How the rest were decided: rule-based first (NCBI synonym of the two names; genus pairs
I listed by hand in `review_disagreements.py`; ANI >= 0.95 across a different lineage; audit-listed top hit; same-genus species at ANI >= 0.99),
then by hand. **26 of the 27 leave-one-out disagreements carry an individual verdict in `data/validation/truth_set_reviews.tsv` (the 27th, a
*Botrytis/Botryotinia* genus pair, is the rule's verdict, which I agree with); in leave-species-out 51 of 91 have an individual verdict and the
others follow the same rules; in leave-genus-out 37 of 410 have one, and the remaining rows are the rules' verdicts, with no individual look
except at the confident ones (status ok), whose list I read.** The 144 leave-genus-out `tool_error` rows are mostly `low_confidence` calls of the
sister family or class (Aspergillaceae/Trichocomaceae, Hypoxylaceae/Xylariaceae, Podosporaceae/Lasiosphaeriaceae, Suillaceae/Auriculariaceae ...);
some of those are taxonomic differences between sources, not errors, and a systematist should look before they are counted against the tool.

**T5. Disagreements (a wrong call at a rank with a known truth) by category**
| dataset | disagreements | tool_error | reference_label_error | taxonomy_synonym | species_complex | low_quality_input | truth_label_error | contamination | status ok (confident) | reviewed by hand |
|---|---|---|---|---|---|---|---|---|---|---|
| loo.known | 27 | 8 | 10 | 2 | 7 | 0 | 0 | 0 | 18 | 26 |
| los.known | 91 | 20 | 41 | 5 | 25 | 0 | 0 | 0 | 70 | 51 |
| log.known | 410 | 144 | 247 | 19 | 0 | 0 | 0 | 0 | 203 | 37 |
| los.novel | 9 | 0 | 9 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| log.novel | 65 | 6 | 59 | 0 | 0 | 0 | 0 | 0 | 0 | 3 |
| indep.known | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |


**loo.known: by rank**
| rank | tool_error | reference_label_error | taxonomy_synonym | species_complex | low_quality_input | total | of which status ok |
|---|---|---|---|---|---|---|---|
| phylum | 1 | 0 | 0 | 0 | 0 | 1 | 0 |
| class | 3 | 1 | 0 | 0 | 0 | 4 | 0 |
| order | 3 | 0 | 0 | 0 | 0 | 3 | 0 |
| family | 1 | 1 | 0 | 0 | 0 | 2 | 1 |
| genus | 0 | 2 | 2 | 0 | 0 | 4 | 4 |
| species | 0 | 6 | 0 | 7 | 0 | 13 | 13 |

**los.known: by rank**
| rank | tool_error | reference_label_error | taxonomy_synonym | species_complex | low_quality_input | total | of which status ok |
|---|---|---|---|---|---|---|---|
| phylum | 2 | 0 | 0 | 0 | 0 | 2 | 0 |
| class | 5 | 5 | 0 | 0 | 0 | 10 | 3 |
| order | 6 | 3 | 0 | 0 | 0 | 9 | 3 |
| family | 6 | 5 | 1 | 0 | 0 | 12 | 6 |
| genus | 1 | 8 | 4 | 0 | 0 | 13 | 13 |
| species | 0 | 20 | 0 | 25 | 0 | 45 | 45 |

**log.known: by rank**
| rank | tool_error | reference_label_error | taxonomy_synonym | species_complex | low_quality_input | total | of which status ok |
|---|---|---|---|---|---|---|---|
| phylum | 10 | 5 | 0 | 0 | 0 | 15 | 0 |
| class | 29 | 64 | 0 | 0 | 0 | 93 | 34 |
| order | 33 | 59 | 0 | 0 | 0 | 92 | 37 |
| family | 63 | 59 | 1 | 0 | 0 | 123 | 45 |
| genus | 9 | 51 | 18 | 0 | 0 | 78 | 78 |
| species | 0 | 9 | 0 | 0 | 0 | 9 | 9 |


**Reference genomes behind the label conflicts.** The 273 label-conflict rows of the known-threshold runs (provisional "sp." names left out) rest on 44 distinct reference genomes. The
largest: `GCA_003069565.1` (labelled *Aspergillus fumigatus*; the best hit of 18 inputs of other classes), `GCA_037044395.1` (*Alternaria destruens*, 7),
`GCA_022813265.1` (*Zalaria obscura*, 6), `GCA_002104975.1` (*Neocallimastix* sp. G1, 4), `GCA_037043785.1` (*Talaromyces purpureogenus*, 4): all
already on the label-audit lists (`reports/label_audit/`). **Fourteen are not on those lists** (sequence similarity only, I have not looked
at the assemblies): `GCA_040938415.1` (*Sarcopodium* sp. QM3-1; best hit of 3 *Fusarium*), `GCA_028028555.1` (*Pseudopestalotiopsis chinensis*
vs *Pestalotiopsis*, ANI 0.991), `GCA_051527695.1` (*Neonectria* sp. MB3469 vs a *Cordyceps*, ANI 0.991), `GCA_008122225.1` (*Ciboria
shiraiana* vs *Sclerotinia sclerotiorum*, ANI 0.997), `GCA_049307545.1` (*Candidozyma auris* vs *Candida parapsilosis*, ANI 0.9985),
`GCA_002116345.1`, `GCA_003325415.1`, `GCA_030782805.1`, `GCA_030247855.2`, `GCA_051529115.1`, `GCA_056321185.1`, `GCA_037044355.1`,
`GCA_937900885.1`, `GCA_041430495.1`. To report upstream under #14.

## Effect of the 43 suspect reference genomes

The production index still contains the 43 genomes of `reports/label_audit/exclude_ids.txt`. Hiding them for every query (`--exclude-ids`):

**T4. Reference genomes on the label-audit list (43) hidden, flag-aware**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out, audit list hidden | 653 | 0.94 / 0.998 [0.991, 1.000] (611) | 0.95 / 0.995 [0.985, 0.998] (603) | 0.94 / 0.995 [0.985, 0.998] (598) | 0.91 / 0.997 [0.987, 0.999] (576) | 0.77 / 0.992 [0.980, 0.997] (503) | 0.39 / 0.949 [0.914, 0.970] (253) |
| leave-species-out, audit list hidden | 653 | 0.88 / 0.997 [0.987, 0.999] (574) | 0.90 / 0.988 [0.975, 0.994] (572) | 0.88 / 0.986 [0.972, 0.993] (562) | 0.83 / 0.979 [0.963, 0.988] (528) | 0.54 / 0.966 [0.941, 0.980] (350) | 0.07 / 0.000 [0.000, 0.080] (44) |
| leave-genus-out, audit list hidden | 653 | 0.66 / 0.965 [0.943, 0.979] (430) | 0.67 / 0.897 [0.865, 0.923] (428) | 0.62 / 0.832 [0.791, 0.865] (392) | 0.50 / 0.722 [0.671, 0.769] (317) | 0.10 / 0.000 [0.000, 0.054] (67) | 0.01 / 0.000 [0.000, 0.354] (7) |


Production index (audit list kept) for comparison:

**T4b. Audit list kept (the production index), flag-aware**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out | 653 | 0.94 / 0.998 [0.991, 1.000] (611) | 0.95 / 0.993 [0.983, 0.997] (603) | 0.94 / 0.995 [0.985, 0.998] (597) | 0.90 / 0.997 [0.987, 0.999] (572) | 0.77 / 0.992 [0.980, 0.997] (500) | 0.39 / 0.948 [0.914, 0.970] (252) |
| leave-species-out | 653 | 0.88 / 0.997 [0.987, 0.999] (575) | 0.90 / 0.983 [0.968, 0.990] (573) | 0.88 / 0.984 [0.970, 0.992] (560) | 0.82 / 0.977 [0.960, 0.987] (523) | 0.53 / 0.962 [0.937, 0.978] (346) | 0.07 / 0.000 [0.000, 0.079] (45) |
| leave-genus-out | 653 | 0.67 / 0.966 [0.944, 0.979] (436) | 0.68 / 0.785 [0.744, 0.821] (433) | 0.62 / 0.766 [0.722, 0.805] (393) | 0.50 / 0.616 [0.561, 0.667] (320) | 0.12 / 0.000 [0.000, 0.047] (78) | 0.01 / 0.000 [0.000, 0.299] (9) |


Leave-one-out and leave-species-out move by at most 0.005. **Leave-genus-out class precision rises from 0.785 to 0.897, order from 0.766 to
0.832, family from 0.616 to 0.722, and the number of leave-genus-out disagreements falls from 410 to 287.** The suspects matter only when the
true relatives are absent and a mislabelled neighbour is the nearest genome. Removing them from the next release is cheap (the build already has
the list); the note of the current release says they were kept deliberately.

## Effect of the decoy panel (#79)

Same taxonomy, index 2026.10.1 (no decoys) against 2026.10.2 (81 decoys), raw tables:

**T3. Index 2026.10.1 (same taxonomy, no decoys), raw**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out without decoys | 653 | 0.98 / 0.978 [0.964, 0.987] (638) | 0.99 / 0.968 [0.951, 0.979] (628) | 0.97 / 0.968 [0.951, 0.979] (619) | 0.92 / 0.985 [0.971, 0.992] (584) | 0.77 / 0.992 [0.980, 0.997] (500) | 0.39 / 0.948 [0.914, 0.970] (252) |
| leave-species-out without decoys | 653 | 0.96 / 0.962 [0.944, 0.974] (626) | 0.97 / 0.942 [0.921, 0.958] (621) | 0.95 / 0.939 [0.917, 0.955] (606) | 0.87 / 0.951 [0.930, 0.966] (551) | 0.53 / 0.962 [0.937, 0.978] (346) | 0.07 / 0.000 [0.000, 0.079] (45) |
| leave-genus-out without decoys | 653 | 0.94 / 0.861 [0.831, 0.886] (612) | 0.95 / 0.686 [0.648, 0.722] (602) | 0.86 / 0.659 [0.619, 0.698] (546) | 0.61 / 0.573 [0.523, 0.621] (386) | 0.12 / 0.000 [0.000, 0.047] (78) | 0.01 / 0.000 [0.000, 0.299] (9) |


**T3b. Index 2026.10.2 (with decoys), raw**
| mode | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| leave-one-out | 653 | 0.98 / 0.950 [0.930, 0.964] (642) | 0.99 / 0.952 [0.933, 0.966] (629) | 0.97 / 0.961 [0.943, 0.974] (618) | 0.92 / 0.971 [0.954, 0.982] (587) | 0.77 / 0.992 [0.980, 0.997] (500) | 0.39 / 0.948 [0.914, 0.970] (252) |
| leave-species-out | 653 | 0.97 / 0.905 [0.880, 0.926] (633) | 0.98 / 0.905 [0.880, 0.926] (622) | 0.94 / 0.923 [0.899, 0.942] (597) | 0.87 / 0.924 [0.899, 0.943] (553) | 0.53 / 0.962 [0.937, 0.978] (346) | 0.07 / 0.000 [0.000, 0.079] (45) |
| leave-genus-out | 653 | 0.95 / 0.679 [0.641, 0.715] (620) | 0.95 / 0.561 [0.521, 0.600] (606) | 0.85 / 0.557 [0.515, 0.599] (540) | 0.68 / 0.454 [0.408, 0.501] (434) | 0.12 / 0.000 [0.000, 0.047] (78) | 0.01 / 0.000 [0.000, 0.299] (9) |


In the raw tables the decoys cost phylum precision: 0.978 to 0.950 in leave-one-out, 0.962 to 0.905 in leave-species-out, 0.861 to 0.679 in
leave-genus-out. The mechanism: a fungus with no close relative has no good hit anywhere, and its nearest genome is by chance a decoy
(*Mus musculus*, *Thalia democratica*, *Lampetra planeri*; ANI 0.75-0.88, a few dozen shared hashes). The phylum vote then reads "Chordata". It
happened to **31 of 653 inputs in leave-one-out, 58 in leave-species-out and 184 (28%) in leave-genus-out, every one `low_confidence` and flagged
`distant_non_fungal_hit`**. Without decoys the phylum call for those inputs was right for 13, 28 and 103 of them and wrong for 13, 22 and 67
(the rest got none): about half to 60% of the calls they would have made were right. Treating them as no call (the flag-aware tables) restores
precision (0.998 / 0.997 / 0.966) but loses the call, hence the lower called fractions. **Suggested change (not made here; `crates/` belongs to
another change):** when the best-hit ANI is below 0.90 and the voted phylum is a decoy phylum, report no phylum (and no lower rank) and keep the
flag; this is what the flag-aware scoring does, and it makes the documented meaning of `distant_non_fungal_hit` true of the call too.

## Outside the index

The 21 RefSeq assemblies whose accession number is not in the index (model organisms such as *Pyricularia oryzae*, *Candida albicans*,
*Trichoderma reesei*: the species is in the index through another assembly). All 21 `ok`; correct at phylum, class, order, family and
genus wherever a truth exists; species called and correct for 17 of 21; the 4 withheld are *Meyerozyma guilliermondii*, *Fusarium verticillioides*,
*Paracoccidioides brasiliensis* and *Trichoderma virens* at ANI 0.9997-0.9999 with `low_agreement` (near-best hits of different species labels:
label inconsistency in the reference, not a failure of the match). Nine are type material.

## Limits

* The truth is the NCBI lineage as of the pinned taxdump; a taxonomy change that NCBI has not made is counted against the tool
  (*Pichiomycetes*, genus splits). Where a name pair is a known synonym or split I categorised it `taxonomy_synonym` (list in
  `review_disagreements.py`, my reading, not an authority).
* The truth set is dominated by well-studied genera; only 21 of 674 are outside the index, so the independent estimate is small.
* "Leave-one-out" hides one accession; identical copies under other accessions remain (14 inputs, no material effect). It is not
  leave-strain-out.
* The leave-genus-out rule hides the genus of the genome as the index holds it; a lineage whose genus is polyphyletic in NCBI keeps relatives.
* Label conflicts are a similarity argument: two genomes at ANI >= 0.95 cannot both be right about a different genus, but which of the two is
  wrong, or whether the assembly is contaminated, was not checked against the assembly. Report, do not delete.
* Disagreement verdicts for leave-genus-out are mostly rule-based (see above); the tables' counts are exact, the verdicts are not.

## Acceptance (#39)

| criterion | status |
|---|---|
| >= 250 type-material genomes scored in all three modes | 258 in leave-one-out, leave-species-out and leave-genus-out |
| every disagreement categorised | yes: 602 rows (the known and novel runs), 4 categories used; 114 of the 528 known-threshold rows carry an individual verdict, all but one of leave-one-out |
| share of BFD-benchmark errors that were label errors | query side 0 (labels identical); reference side 37% / 45% / 60% of the disagreements (above) |
| accuracy on the trusted set next to the full-BFD benchmark | above |
