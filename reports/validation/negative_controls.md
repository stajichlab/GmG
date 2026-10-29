# Negative controls: real non-fungal genomes (#26)

`scripts/hpcc/16_negative_controls.sbatch`: seven NCBI RefSeq genomes classified against the BFD index (balanced default thresholds, `gmg run`). Expected result: no fungal call.

| genome | group | size | status | deepest rank called | top hit (species) | ANI | shared hashes / query hashes | agreement on phylum |
|---|---|---|---|---|---|---|---|---|
| *Escherichia coli* K-12 (GCF_000005845.2) | bacterium | 4.6 Mb | `no_match` | none | *Seiridium cardinale* | 0.775 | 2 / 424 | 1.00 |
| *Bacillus subtilis* 168 (GCF_000009045.1) | bacterium | 4.2 Mb | `no_match` | none | *Gigaspora rosea* | 0.793 | 3 / 391 | 1.00 |
| phiX174 (GCF_000819615.1) | virus | 5.4 kb | `no_match` | none | none | n/a | no shared hash | n/a |
| *Arabidopsis thaliana* (GCF_000001735.4) | plant | 120 Mb | `low_confidence` | order | *Gigaspora margarita* | 0.799 | 97 / 10,704 | 0.68 |
| *Chlamydomonas reinhardtii* (GCF_000002595.2) | green alga | 111 Mb | `low_confidence` | class | *Purpureocillium* sp. | 0.759 | 28 / 9,270 | 1.00 |
| *Caenorhabditis elegans* (GCF_000002985.6) | nematode | 100 Mb | `low_confidence` | class | *Gigaspora margarita* | 0.807 | 99 / 9,037 | 0.52 |
| *Plasmodium falciparum* 3D7 (GCF_000002765.5) | apicomplexan | 23 Mb | `low_confidence` | family | *Neocallimastix lanati* | 0.899 | 203 / 1,902 | 1.00 |

Result: **bacteria and the virus pass** (`no_match`, with the `few_hashes` situation of tiny inputs for phiX). **Every eukaryote fails**: a fungal phylum is named (Glomeromycota for the plant and the nematode, Ascomycota for the alga, Neocallimastigomycota for the malaria parasite) and the class, order or family is called for some. No genus or species is called, because the genus cutoff (ANI 0.90) is out of reach, but the phylum threshold (ANI 0.7109, i.e. 0.08% of the query's hashes shared) is far below what a 100 Mb genome shares with *any* fungal genome by chance and conserved sequence (organelle, rDNA, repeats, AT-rich composition: the *Plasmodium* genome shares 10.7% of its hashes with an AT-rich *Neocallimastix*).

Why: the calibration (`reports/calibration_known.md`, `novel`) was fitted on queries that are fungi; nothing in the reference tells the caller what a non-fungus looks like. The agreement of the near-best hits and the shared-hash count do not separate these cases from distant fungi (a genuine Microsporidium vs other fungi sits at the same ANI 0.75-0.80). The human genome (3 GB) was not run.

Options (a decision for the maintainers, see the issue): add decoy genomes from non-fungal lineages (plants, animals, protists, bacteria, archaea, viruses, human) to the reference as an explicit "not a fungus" group, so a non-fungal query's best hits are decoys and the answer is a stated `not_fungal` rather than a wrong fungal phylum; and/or calibrate a floor on shared hashes relative to genome size for the phylum call. Until then the documentation says: GmG assumes the input is a fungus (docs/tutorial.md, section 4).
