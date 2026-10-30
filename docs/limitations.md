# What GmG can and cannot say

GmG compares the k-mer content of an input with a reference collection and reports how similar the input is to the closest reference genomes. Everything below is a measurement on the reference and datasets listed in [reports/validation/README.md](../reports/validation/README.md) (BFD index 2026.10.x, 23,669 genomes, scaled = 10000, default `known` thresholds). Numbers will change with the reference and are re-measured for each index release.

## What a result means

* A call at a rank (phylum to species) means the input's nearest reference genomes agree on that taxon and are similar enough (ANI cutoffs per rank). It is a similarity statement, not a taxonomic opinion. A name is only as good as the reference label (some BFD/NCBI labels are wrong, see `reports/label_audit/`).
* `status = ok` means a call down to genus or species was made; `low_confidence` means only a shallow call or a `tentative_genus`; `no_match` means nothing in the reference is close enough. A `tentative_genus` is never a call.
* Species names follow the reference label. Species complexes and recently split species are not resolved by a sketch: two close species of one genus (for example *Coccidioides immitis* and *C. posadasii*, ANI about 0.976-0.980) are told apart only when the species is in the reference.
* Hits at 100% ANI (`identical_to_reference`) mean the input, or a near-identical copy, is in the reference.

## Measured accuracy (calls made / precision)

| situation | measure | result |
|---|---|---|
| RefSeq truth set, a new strain of a known species (leave-one-out, 653 genomes, 258 type material; index 2026.10.2) | precision of calls made, phylum / class / order / family / genus | 0.998 / 0.993 / 0.995 / 0.997 / 0.992 (species 0.948, called for 39%; 0.992 and 60% when the species has another genome in the reference) |
| same, a new species of a known genus (leave-species-out) | phylum to genus | 0.997 / 0.983 / 0.984 / 0.977 / 0.962; a species call here is necessarily wrong and happens for 7% (congeners at ANI >= 0.99, label conflicts) |
| same, a new genus of a known family (leave-genus-out) | phylum / class / order / family | 0.966 / 0.785 / 0.766 / 0.616 (called 67% / 68% / 62% / 50%); 12% get a genus call (label conflicts, synonym genera, 9 plain overcalls). The 57 ambiguous-label genomes left out from 2026.10.3 raise class to about 0.90 (measured with the 43 audit genomes hidden: class 0.897, order 0.832, family 0.722) |
| RefSeq genomes outside the reference (21) | phylum to genus | 21 / 21 correct at each rank; species called and correct for 17 |
| Independent assemblies (1,257 not in the reference, no quality filter), balanced default | precision of calls made | phylum 0.990, class 0.991, order 0.986, family 0.983, genus 0.988 |
| same set | share with a call | phylum 0.946; genus 0.839; inputs with only a `tentative_genus` 78 |
| Species (52 genomes with a species name, all in the reference) | called / precision | 75% / 1.000 (too few to generalise) |
| Genus absent from the reference (22 genomes) | wrongly named at genus | 3 of 22 (leave-genus-out in BFD: 18%) |
| Species absent from the reference | wrongly named at species | 11.7% of queries |
| Real reads, *Coccidioides* (558 sets, 100k to 2M reads) | genus / species called, precision | genus ≥ 99.6% called, 0.991-0.993; species 95-96% called, 0.994-0.996 |
| Real RNA-seq (345 sets) | genus / species precision | genus 0.997 at every depth; species 0.966-0.968, called 63-68% |
| Degraded genomes (24 hosts) | completeness 10-75%, 5 kb contigs | calls unchanged down to 10% (phylum right for 23/24) |
| Contamination from another phylum | phylum correct | 21/24 at 5%, 19/24 at 20%, 17/24 at 50% foreign DNA (another fungus) |
| Bacterial (to 15%) or human (to 30%) DNA in an assembly (24 + 44 hosts) | genus and phylum call unchanged; species | genus/phylum unchanged for 255/256 and 474/474; species kept 31/32 at 5%, 29-30/32 at 10%, withheld (never wrong) from 15% bacteria / 20% human |
| Negative controls | random DNA, repeat, bacteria, phage | `no_match` |
| Non-fungal eukaryotes | see below | flagged `not_fungal` when a relative (ANI >= 0.90) is among the decoys; otherwise `distant_non_fungal_hit` |

Precision is the share of calls made that are correct; "called" is the share of inputs that get a call. Confidence intervals and per-clade tables are in the validation reports.

## Where it is weaker

* **Clades that are thin in the reference** (small phyla such as Neocallimastigomycota, Basidiomycota genus precision 0.942 before the balanced default): fewer calls, and calls are more often shallow. Clade-aware calibration: #58.
* **Genera and species that are not in the reference** are sometimes given the nearest name. The `low_confidence` status and `tentative_genus` exist for this; do not treat an ANI-0.9 genus call as a new-taxon discovery or a negative result.
* **Small, fragmented or partial inputs** (N50 under 10 kb, fewer than 100 hashes) often get no call; this is correct but says nothing about the organism. QC flags `fragmented`, `few_hashes`, `short_genome`, `long_genome`, `high_N` say why.
* **Mixtures.** GmG returns one lineage. A contaminant that has close relatives in the reference can win the call; a mixture flag was tried and rejected (6% false positives at ~50% sensitivity). GmG does not separate the parts of a metagenome (#31).
* **Not a fungus.** The input is assumed to be a fungus. A decoy panel of non-fungal genomes (plants, animals, protists) lets non-fungal inputs that have a close decoy relative (best-hit ANI >= 0.90) say so (`not_fungal` flag, #79). Inputs far from every reference genome, fungal or not, usually land on whichever decoy is nearest by chance (ANI 0.75-0.88); they get `distant_non_fungal_hit`; the call is then made without the decoys, so it is a low-confidence call (as without decoys) that says nothing about whether they are fungi (about 2% of divergent real fungi in the post-freeze set look like this). Bacteria and viruses return `no_match`.
* **Reads** need enough depth for the species to be seen at the species cutoff (0.99 for reads, 0.9938 for assemblies): at about 0.5x genome coverage genus calls are reliable, species calls less so. Sequencing errors lower ANI. RNA-seq covers only the transcribed part of the genome, so it is less sensitive (#33). Long reads and very high error rates are not validated yet (#44).
* **Only fungi are in the reference.** Mitochondrion-only or plasmid-only inputs are not tested.

## What GmG is not

It is not a substitute for taxonomic curation, a phylogenetic placement, an ANI-based species delimitation, or a contamination or completeness assessment (use BUSCO, CheckM2, FCS-GX or similar). Use the call as a fast first answer and a prompt for a check, and report the index version (`gmg info`; recorded in `run.json`) with it.

## Reference genomes left out for ambiguous labels

57 genomes whose label conflicts with their close relatives (43 from the label audit, `reports/label_audit/`, and 14 more found by the truth-set run, #39) are
not in the reference from release 2026.10.3 on; the list is `reports/label_audit/reference_exclusions.txt`. Hiding the 43 audit genomes raised
leave-genus-out class precision from 0.785 to 0.897 (`reports/validation/truth_set.md`). A species whose only genome was excluded loses its species call.

## Contamination

Genus and phylum calls are unchanged by up to 15% bacterial or 30% human DNA in an assembly (`reports/validation/contamination.md`). The species call needs about 90% of the input to be fungal DNA of the species: ANI falls as (1 - f)^(1/21), and the species cutoff (0.9938) is crossed at 10-15% foreign DNA, after which the species is withheld (never wrong) and the genus is still called.
