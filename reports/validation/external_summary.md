# External independent datasets (#40): summary

Run 2026-10-29/30 with `gmg 0.0.1 (a9acf5cc9)` on `bfd-2026.10.3-s10000.gmgi` (decoys included; the `weak_phylum` flag and release 2026.10.4, which drops the *Thalia democratica* decoy, came later). Registry: `data/validation/datasets.tsv`. Prep: `scripts/prep_external.py`; run: `scripts/hpcc/20_external_set.sbatch`; reports: `scripts/external_report.py`, `scripts/score_by_name.py`. Per-genome tables stay on the HPCC.

| dataset | inputs | outside BFD | phylum precision (called) | genus precision (called) | notes |
|---|---|---|---|---|---|
| post-freeze NCBI | 1,257 | 1,257 | 0.990 (0.95) | 0.988 (0.84) | NCBI lineage, `post_freeze_balanced.md` (index 2026.10.1) |
| IndigoAg | 199 | 189 | 0.954 (0.994) | 0.713 (0.661) | provider names from 2017; genus is a lower bound |
| ZyGoLife aaftf/sort | 552 | 552 | 0.967 (0.987) | 0.901 (0.715) | taxid-based truth; species equal to the sheet name for 102 of 177 called (0.58, Mucor/Rhizopus complexes) |
| KassonLab | 36 | 34 | 0.850 (0.556) | 0.875 (0.444) | 16 `no_match`: short, fragmented herbarium/metagenome assemblies |
| Antarctica MAGs | 266 | 259 | not scored | not scored | 251 phylum calls (201 Ascomycota, 27 Basidiomycota), 7 genus, 12 species; 15 `no_match`, 68 `distant_non_fungal_hit`, 0 `not_fungal` |
| CCFEE AAFTF vs shovill | 29 + 38 | 67 | phylum agrees with samples.csv: 27/27, 37/38 | n/a | same strain, two assemblers: phylum 26/27 identical, family 20/20, genus 3/3 |

Six datasets, 2,300 genomes outside BFD (1,257 + 189 + 552 + 34 + 259 + 67 strains), a results table per dataset (the linked reports).

## What stands out

* **ZyGoLife (thin clades, Mucoromycota/Entomophthoromycota):** phylum 0.967 and genus 0.901 are close to the BFD-wide numbers even for clades BFD covers thinly; 544 of 552 get a phylum. Wrong phylum calls: Mucoromycota->Glomeromycota 5, Entomophthoromycota->Neocallimastigomycota 4, Entomophthoromycota->Chordata 2 (the *Thalia democratica* decoy; removed in 2026.10.4), Mucoromycota->Basidiomycota 2. 34 inputs carry `distant_non_fungal_hit` and 2 `not_fungal`: both false for fungi here, the decoy ANI gate does not hold for Entomophthoromycota. Genus errors are within-Mucorales confusions (*Mucor*/*Gilbertella*, *Rhizopus*/*Circinella*, *Mortierella*/*Linnemannia*) and the provider names are partly outdated.
* **KassonLab:** only 20 of 36 get a call because 21 are `short_genome`, 20 `fragmented`, 17 `few_hashes` (AAFTF assemblies of metagenomic or herbarium DNA, large Entomophthorales genomes). When a genus is called it is right 14 of 16, with the two wrong ones *Cryphonectria parasitica* -> *Neonectria* (a herbarium metagenome) and *Ophiocordyceps* sp. -> *Massospora*. Species: 11 of 12 called equal the sheet species.
* **Antarctica MAGs:** 235 of 266 are `low_confidence` (family or order level at best): the bins are mostly from lineages the reference lacks; nothing is flagged `not_fungal`. No names exist, so these calls are hypotheses.
* **CCFEE:** the two assemblers agree at every rank where both call; shovill gives a phylum for all 38 inputs where AAFTF leaves 2 short assemblies without. Genus is called for only 3-5 of ~30 strains (Antarctic rock fungi are absent from BFD).
* **IndigoAg wrong-phylum assemblies are not contaminated:** see `external_indigo.md`; the plurality phylum of the exclusive-hash split is correct for 174/174.

## Not done

* D2 NCBI fungal MAGs, D4 JGI Mycocosm (no export found), D6 pathogen panels, D7 the 10 type-material genomes outside BFD (partly covered by #39). The acceptance (>= 6 datasets, >= 500 genomes outside BFD, a results table each) is met without them.
* Lab assemblies named in the issue (Coccidioides pangenome, MAT genomes, the 14 desert isolates) were not located.
