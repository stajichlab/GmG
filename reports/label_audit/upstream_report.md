# BFD and lab label hand-off report (GmG #14)

**Status: internal hand-off document, written 2026-10-30. Nothing has been sent to anyone.** No recipient outside the lab exists yet; NCBI is a possible later recipient, so the
table is laid out so that a copy can be cut down to *accession, current organism name, suggested taxon, evidence* (columns of `upstream_report.tsv`). No external issue was opened and no email was written.

**What the suspicions are based on.** Every suspicion here comes from **sequence similarity only** (k-mer sketches, containment ANI, 95%-ANI clusters, and the labels of the
neighbouring genomes). None has been checked against the biological material, the submitter, or a curated taxonomy, and similarity cannot say which of two conflicting
genomes carries the wrong name. **The BFD is the lab's own reference** (about 23.7k fungal genomes with the lineage derived from the NCBI taxid of each assembly record at build time),
not an authority: "current organism" below is the BFD species label, and the NCBI record may have been corrected since. Check the current name on NCBI before anything is forwarded.
"ANI" is `containment(query)^(1/21)` of the scaled=10000 sketch, so 1.0000 means every sketch hash of the query is in the neighbour; it is not a classical ANI.

## Counts

<!-- COUNTS -->
| group | A | B | C | total |
|---|---|---|---|---|
| BFD genomes left out: label audit (43) | 20 | 6 | 17 | 43 |
| BFD genomes kept in 2026.10.3: conflict with another kept genome (the 9 A and 6 B left out from 2026.10.5) | 9 | 6 | 70 | 85 |
| BFD genomes left out: truth-set run (14) | 9 | 3 | 2 | 14 |
| Lab or provider labels of read sets / assemblies (not NCBI accessions) | 5 | 9 | 9 | 23 |
| SRA run (ERR14883901, assembly of the same BioSample) | 0 | 1 | 0 | 1 |
| BFD genomes kept: size outliers only | 0 | 0 | 69 | 69 |
| Decoy assembly GCF_965202585.1 | 0 | 0 | 1 | 1 |
| **all rows** | **43** | **25** | **168** | **236** |
<!-- /COUNTS -->

* **A, strong**: at least two independent lines of evidence and the genome is the odd one out (details under "Tiers"). Suggested action is relabel or flag.
* **B, moderate**: one line at ANI >= 0.98, or two lines that do not make it to A, or a pair of genomes where the data do not say which side is wrong. Suggested action is "check".
* **C, weak or ambiguous**: uninformative `sp.` labels, ANI < 0.98 from a single line, order-level conflicts, size outliers only, metagenome inputs, taxonomy updates. Suggested action is "keep".

Where the genome sits matters more than the tier for a copy to NCBI. Of the 236 rows, **57 genomes are left out of the GmG reference from index 2026.10.3** (`reference_exclusions.txt`, lines 1-57:
the 43 of the label audit plus 14 found by the truth-set run, #39/#90) and **15 more from index 2026.10.5** (lines 58-72: the 9 tier A and 6 tier B BFD genomes that were still in the reference when this report was written;
the maintainer decided on 2026-10-30 to leave them out, as was done for the 57; column `in_reference_2026.10.5` of `upstream_report.tsv`). **138 BFD genome rows remain in the reference** (70 tier C label conflicts with another kept genome, 68 size outliers; the 69th, GCA_003314755.1, is one of the 15;
the SRA genome pair is a further 2) and 24 rows are outside BFD (23 lab or provider label rows and one decoy). Deliberately **not** left out: the tier C rows (weak or ambiguous evidence)
and the assembly GCA_982397225.1 of the SRA pair ERR14883901 (which of the read set and the assembly carries the wrong name is undetermined). The row counts and tiers below are those of the report as written against 2026.10.3; the tables list the 15 as still in the reference of that release.

Highlights, all to be read with the caveat above:

* **29 of the 57 left-out genomes are tier A, 9 tier B and 19 tier C.** The A rows are the clear cases: e.g. GCA_049307695.1 "Emarellia grisea" is the nearest genome of 416 genomes labelled Eurotiomycetes and sits in a 95%-ANI cluster of 415 others, 414 of them *Aspergillus fumigatus*; GCA_037044395.1 "Alternaria destruens" sits in a cluster with 81 other *Aspergillus* genomes (82 members besides itself).
* **19 of the 57 are probably not wrong**, they are tier C: twelve carry only an uninformative label (mostly "Saccharomycetales sp." whose neighbours are in the new class Pichiomycetes), five are specific names probably better labelled than their neighbour (the RefSeq *Brettanomyces bruxellensis* GCF_011074885.1, *Candida tropicalis* GCA_057383725.1, *Alternaria* sp. SPS-2, *Tetracladium maxilliforme*, *Brettanomyces* sp. HC-2020a) and two are single weak lines (ANI < 0.98). They were left out as ambiguous, which is conservative; this report suggests re-admitting nothing, but the lab may want to look.
* **Surprise: 15 genomes that were still in the 2026.10.3 reference look as wrong as the left-out ones (left out from 2026.10.5)** (tier A or B, flag `kept_in_reference`). The audit defined hubs only at phylum/class level; at order level the same pattern exists, e.g. GCA_011057605.1 "Pseudocercospora fijiensis" (nearest genome of 19 Botryosphaeriales genomes; its 43 cluster mates are 41 Botryosphaeriales, mostly *Lasiodiplodia*), GCA_023087585.1 "Bipolaris maydis" and GCA_023079265.1 "Alternaria brassicicola" in the same cluster, GCA_037044435.1 "Alternaria alstroemeriae" (ANI 1.0 to an *Aspergillus costaricensis*, in the accession range GCA_037044xxx that also holds five left-out genomes). These are not in `reference_exclusions.txt`.
* **ERR14883901** (BioSample SAMEA118140668, registered as *Coccidioides immitis*) and its assembly GCA_982397225.1: best hit *Ophidiomyces ophidiicola* GCA_982397405.1 at ANI 0.947 (reads, own genome hidden) and 0.950 (assembly), not any of 5 other *C. immitis* in the reference. The assembly is in the reference (kept: the wrong side is undetermined). Tier B (same BioSample, so the lines are not fully independent).
* **Lab and provider labels** (not NCBI accessions): the strongest are four strains whose lab name disagrees with both the reads and the assembly at ANI 0.986-1.0 (*Absidia* sp. A-19777 is *Lichtheimia ramosa*, *Gongronella* sp. A-20010 is *Halteromyces radiatus*, *Rhizomucor pusillus* 2543 is *Phycomyces nitens*, *Chaetocladium brefeldii* RSA 1136 is *Cunninghamella blakesleeana*), plus *Nowakowskiella* sp. JEL0774, whose reads match a reference genome that carries the same strain designation under the name *Gonapodya* (ANI 0.990).
* **The *Thalia democratica* decoy assembly GCF_965202585.1** shares 11-15% of its k-mers with *Conidiobolus*-type fungi, closer than any fungal reference genome; possible fungal contamination of that NCBI assembly (not examined; tier C).

## Methods (short)

The 43 + 23 + 69 audit lists (`find_label_suspects.py`, `analyze_clusters.py`; leave-one-out of all 23.7k BFD genomes at scaled=10000, k=21): a genome is a suspect at phylum, class or order when the vote of its
nearest genomes (within 10% of the best shared count, consistency >= 0.8) names a different taxon and the best-hit ANI is >= 0.95; a hub is the nearest genome of >= 3 conflicting phylum/class queries; the 43
excluded ids are the phylum/class suspects plus the hubs; the 69 size outliers have a sketch more than 1.6x the median of their species (species with >= 5 genomes) and were kept. For this report each of the 57 left-out genomes, and every
kept genome that conflicts with another kept genome (`loo_suspects.tsv`, `log_conflicts.tsv`, partner not left out; a genome whose conflicting neighbour is one of the 57 is a victim of that neighbour and not listed), was re-examined with
(1) its 95%-ANI cluster (`data/clusters/bfd_ani0.95.tsv.gz`): does any other member share its class (genus for the 14 truth-set exclusions), is it the odd one out; (2) the hub count; (3) RefSeq genomes of the
truth set (`data/validation/truth_set.list.tsv`) that had it as their nearest neighbour in the first truth-set run (index 2026.10.2, `disagreements.tsv` at commit 374dfec; the re-run on 2026.10.3 no longer has the 14 later exclusions as neighbours); (4) the neighbour's own
cluster. Lab and provider labels come from GmG runs of ZyGoLife (552 assemblies and 119 read sets at 1.5M pairs), IndigoAg (199), KassonLab (36) and the SRA reads (100 runs); a label is listed only where the call is made (`ok`) at ANI >= 0.95 and differs at order or higher (assemblies), or reads and assembly agree against the label. Per-genome files of those runs stay on the HPCC (`GmG_runs/external*`, `reads/titration`, `reads_sra`); the rows were read from them, read-only. The table is regenerated by `scripts/make_upstream_report.py`.

## Tiers, exactly

* **A**: the genome has a specific name (genus given), is the odd one out of its 95%-ANI cluster (>= 2 other members, a single coherent other class/order/genus among them, none or <= 5% sharing its own) or its nearest neighbour is supported by its own cluster, **and** at least one more line: it is a hub (>= 3 other genomes point at it), RefSeq genomes of another lineage point at it, or the ANI is >= 0.98. For the lab rows: reads and assembly of the same strain agree against the label at ANI >= 0.98, or reads match a genome with the same strain designation.
* **B**: ANI >= 0.98 (or two lines) but the side is not determined (the genome and its partner are alone or tied in their clusters), or a phylum-level conflict at ANI >= 0.95; for lab rows one assembly or read set alone at ANI >= 0.98 / a class change, or reads and assembly agree below 0.98.
* **C**: the label is only a higher rank or `sp.` ("uninformative"), the label is supported by its cluster, ANI < 0.98 from a single line, a metagenome or contaminated-looking input, a size outlier with no other finding, or a taxonomy update.

## Tables

Columns: id; accession (NCBI accession, or the lab id for lab rows); current organism (BFD species label, or the lab name); GmG best hit (accession and taxon of the nearest genome); ANI; size (sketch size over the species median, only for the
size outliers); flags; evidence (LOO = leave-one-out conflict, HUB, REFSEQ = RefSeq genome of another lineage, CLUSTER, SIZE, READS, ASM, SRA, STRAIN); suggested taxon. The full evidence sentence, the suggested action, the lineage and the sources are in `upstream_report.tsv`.

<!-- TABLES -->
### Tier A

#### BFD genomes (57 left out of the reference from 2026.10.3 and 15 more from 2026.10.5; the others kept, flag kept_in_reference or size_outlier) (38)

| id | accession | current organism (BFD/NCBI) | GmG best hit | ANI | size | flags | evidence | suggested |
|---|---|---|---|---|---|---|---|---|
| A001 | GCA_002233575.1 | Geotrichum fermentans | GCA_003116895.1 Cutaneotrichosporon dermatis | 1.0000 |  | audit43;cluster95_others=0;neighbour_in_reference | LOO | Cutaneotrichosporon dermatis |
| A002 | GCA_003069565.1 | Aspergillus fumigatus | GCA_024606205.1 Alternaria sp. SPS-2 | 1.0000 |  | audit43;hub_n=210;cluster95_others=209 | LOO+HUB+REFSEQ+CLUSTER | Alternaria alternata (145 of 209) |
| A003 | GCA_022813755.1 | Cladosporium austrohemisphaericum | GCA_022813565.1 Trichoderma citrinoviride | 1.0000 |  | audit43;hub_n=6;cluster95_others=5;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Trichoderma citrinoviride (5 of 5) |
| A004 | GCA_035928285.1 | Vermiconidia calcicola | GCA_035928545.1 Exophiala xenobiotica | 1.0000 |  | audit43;hub_n=34;cluster95_others=34;neighbour_in_reference | LOO+HUB+CLUSTER | Exophiala xenobiotica (33 of 34) |
| A005 | GCA_037044025.1 | Cladosporium antarcticum | GCA_037043915.1 Fusarium verticillioides | 1.0000 |  | audit43;hub_n=77;cluster95_others=70;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Fusarium verticillioides (54 of 70) |
| A006 | GCA_037044395.1 | Alternaria destruens | GCF_003184535.1 Aspergillus eucalypticola CBS 122712 | 1.0000 |  | audit43;hub_n=81;cluster95_others=82 | HUB+REFSEQ+CLUSTER | Aspergillus (81 of 82; e.g. Aspergillus niger) |
| A007 | GCA_046252445.1 | Candidozyma auris | GCA_003116565.1 Aspergillus fumigatus | 1.0000 |  | audit43;cluster95_others=0;neighbour_in_reference | LOO | Aspergillus fumigatus |
| A008 | GCA_001659765.1 | Wickerhamiella versatilis | GCA_979262865.1 Zygosaccharomyces rouxii | 0.9997 |  | audit43;hub_n=43;cluster95_others=23;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Zygosaccharomyces rouxii (21 of 23) |
| A009 | GCA_022813265.1 | Zalaria obscura | GCA_022814395.1 Aspergillus calidoustus | 0.9992 |  | audit43;hub_n=10;cluster95_others=10;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Aspergillus calidoustus (8 of 10) |
| A010 | GCA_049307695.1 | Emarellia grisea | GCA_051225625.1 Aspergillus fumigatus | 0.9991 |  | audit43;hub_n=416;cluster95_others=415;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Aspergillus fumigatus (414 of 415) |
| A011 | GCA_039700845.1 | Alternaria alternata | GCA_033060795.1 Hypocreales sp. ML197-S3 | 0.9982 |  | audit43;cluster95_others=2;neighbour_in_reference | LOO+CLUSTER | Lecanicillium (1 of 2; e.g. Lecanicillium sp. S3808) |
| A012 | GCA_001074275.1 | Naumovozyma dairenensis | GCA_029931525.1 Candida albicans | 0.9975 |  | audit43;hub_n=135;cluster95_others=133;neighbour_in_reference | LOO+HUB+CLUSTER | Candida albicans (132 of 133) |
| A013 | GCA_037043785.1 | Talaromyces purpureogenus | GCA_037893235.1 Trichoderma afroharzianum | 0.9963 |  | audit43;hub_n=18;cluster95_others=18;neighbour_in_reference | LOO+HUB+CLUSTER | Trichoderma afroharzianum (9 of 18) |
| A014 | GCA_025802115.2 | Bertia moriformis | GCA_022813405.1 Aspergillus pseudoglaucus | 0.9956 |  | audit43;hub_n=3;cluster95_others=3;neighbour_in_reference | LOO+HUB+CLUSTER | Aspergillus pseudoglaucus (2 of 3) |
| A015 | GCA_037044075.1 | Phialemonium inflatum | GCA_019022425.1 Talaromyces purpureogenus | 0.9950 |  | audit43;hub_n=8;cluster95_others=8;neighbour_in_reference | LOO+HUB+CLUSTER | Talaromyces purpureogenus (5 of 8) |
| A016 | GCA_037044295.1 | Aspergillus subversicolor | GCA_018350245.1 Fusarium proliferatum | 0.9932 |  | audit43;hub_n=92;cluster95_others=83;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Fusarium (83 of 83; e.g. Fusarium proliferatum) |
| A017 | GCA_019924265.1 | Aureobasidium melanogenum | GCA_043234805.1 Exophiala dermatitidis | 0.9901 |  | audit43;hub_n=38;cluster95_others=36;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Exophiala dermatitidis (35 of 36) |
| A018 | GCA_055690705.1 | Omphalotus olearius | GCA_023619845.1 Aspergillus neotritici | 0.9894 |  | audit43;cluster95_others=0;neighbour_in_reference | LOO | Aspergillus neotritici |
| A019 | GCA_057978835.1 | Cladosporium sp. 450A | GCA_027574645.1 Fusarium solani | 0.9874 |  | audit43;hub_n=4;cluster95_others=2;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Fusarium vanettenii (1 of 2) |
| A020 | GCA_004959885.1 | Penicillium sp. W3 MMC-2018 | GCA_000587855.1 Mucor irregularis | 0.9862 |  | audit43;hub_n=4;cluster95_others=2;neighbour_in_reference | LOO+HUB+CLUSTER | Mucor irregularis (2 of 2) |
| A021 | GCA_037892935.1 | Phlebopus portentosus | GCA_040938115.1 Rigidonotus glomeratus | 1.0000 |  | kept_in_reference;cluster95_others=5;neighbour_in_reference | LOO+CLUSTER | Sanghuangporus sanghuang (3 of 5) |
| A022 | GCA_037044435.1 | Alternaria alstroemeriae | GCA_037044335.1 Aspergillus costaricensis | 1.0000 |  | kept_in_reference;cluster95_others=82;neighbour_in_reference | LOO+CLUSTER | Aspergillus (81 of 82; e.g. Aspergillus niger) |
| A023 | GCA_023087585.1 | Bipolaris maydis | GCA_023087545.1 Lasiodiplodia pseudotheobromae | 0.9995 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Lasiodiplodia (36 of 43; e.g. Lasiodiplodia theobromae) |
| A024 | GCA_049307545.1 | Candidozyma auris | GCF_000182765.1 Lodderomyces parapsilosis | 0.9985 |  | truthset14;cluster95_others=97 | REFSEQ+CLUSTER | Lodderomyces parapsilosis (97 of 97) |
| A025 | GCA_008122225.1 | Ciboria shiraiana | GCF_000146945.2 Sclerotinia sclerotiorum 1980 UF-70 | 0.9969 |  | truthset14;cluster95_others=36 | REFSEQ+CLUSTER | Sclerotinia sclerotiorum (36 of 36) |
| A026 | GCA_937872655.1 | Leucosporidium creatinivorum | GCA_054785315.1 Glaciozyma watsonii | 0.9960 |  | kept_in_reference;hub_n=3;cluster95_others=3;neighbour_in_reference | LOO+HUB+CLUSTER | Glaciozyma watsonii (3 of 3) |
| A027 | GCA_037044355.1 | Aspergillus foetidus | GCF_028827155.1 Penicillium citrinum | 0.9950 |  | truthset14;cluster95_others=35 | REFSEQ+CLUSTER | Penicillium citrinum (30 of 35) |
| A028 | GCA_011057605.1 | Pseudocercospora fijiensis | GCA_018153875.1 Lasiodiplodia theobromae | 0.9938 |  | kept_in_reference;hub_n=19;cluster95_others=43;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Lasiodiplodia (36 of 43; e.g. Lasiodiplodia theobromae) |
| A029 | GCA_023079265.1 | Alternaria brassicicola | GCA_029169085.1 Lasiodiplodia pseudotheobromae | 0.9921 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Lasiodiplodia (36 of 43; e.g. Lasiodiplodia theobromae) |
| A030 | GCA_030247855.2 | Knufia peltigerae | GCF_000835515.1 Exophiala oligosperma | 0.9913 |  | truthset14;cluster95_others=4 | REFSEQ+CLUSTER | Exophiala oligosperma (3 of 4) |
| A031 | GCA_022813475.1 | Pseudopyrenochaeta terrestris | GCA_022814635.1 Cladosporium colombiae | 0.9907 |  | kept_in_reference;hub_n=3;cluster95_others=3;neighbour_in_reference | LOO+HUB+CLUSTER | Cladosporium colombiae (2 of 3) |
| A032 | GCA_051529115.1 | Neohyperdermium piperis | GCF_028009165.1 Akanthomyces muscarius | 0.9846 |  | truthset14;cluster95_others=16 | REFSEQ+CLUSTER | Akanthomyces lecanii (11 of 16) |
| A033 | GCA_017499575.1 | Sanguinoderma rugosum | GCA_024195665.1 Pleurotus ostreatus | 0.9820 |  | kept_in_reference;hub_n=9;cluster95_others=9;neighbour_in_reference | LOO+HUB+CLUSTER | Pleurotus ostreatus (6 of 9) |
| A034 | GCA_002116345.1 | Mycosphaerella sp. PB-2012b | GCF_000320565.1 Sphaerulina musiva SO2202 | 0.9810 |  | truthset14;cluster95_others=18 | REFSEQ+CLUSTER | Sphaerulina musiva (18 of 18) |
| A035 | GCA_040938415.1 | Sarcopodium sp. QM3-1 | GCF_000149955.1 Fusarium oxysporum f. sp. lycopersici 4287 | 0.9785 |  | truthset14;cluster95_others=809 | REFSEQ+CLUSTER | Fusarium oxysporum (782 of 809) |
| A036 | GCA_030782805.1 | Aschersonia sp. MBC 622 | GCF_000280675.1 Beauveria bassiana ARSEF 2860 | 0.9779 |  | truthset14;cluster95_others=227 | REFSEQ+CLUSTER | Beauveria bassiana (221 of 227) |
| A037 | GCA_046127395.1 | Dichomera gemmicola | GCA_041146855.1 Gemmamyces piceae | 0.9645 |  | kept_in_reference;hub_n=5;cluster95_others=6;neighbour_in_reference | LOO+HUB+CLUSTER | Gemmamyces piceae (6 of 6) |
| A038 | GCA_003325415.1 | Rhizopus stolonifer | GCF_025201335.1 Gilbertella persicaria | 0.9609 |  | truthset14;cluster95_others=3 | REFSEQ+CLUSTER | Gilbertella persicaria (3 of 3) |

#### Lab or provider labels of read sets and assemblies, the SRA run, the decoy (5)

| id | accession | current organism (BFD/NCBI) | GmG best hit | ANI | size | flags | evidence | suggested |
|---|---|---|---|---|---|---|---|---|
| A039 | (lab id) NRRL A-19777 | Absidia sp. NRRL A-19777 | GCA_000945115.1 Lichtheimia ramosa | 0.9991 reads; 1.0000 assembly |  | ok (both) | READS+ASM | Lichtheimia ramosa |
| A040 | (lab id) NRRL A-20010 | Gongronella sp. NRRL A-20010 | GCF_025201355.1 Halteromyces radiatus | 0.9963 reads; 0.9997 assembly |  | ok (both) | READS+ASM | Halteromyces radiatus |
| A041 | (lab id) JEL0774 | Nowakowskiella sp. JEL0774 | GCA_027604505.1 Gonapodya sp. JEL0774 (Monoblepharomycota) | 0.9903 |  | ok; 3214 shared hashes | READS+STRAIN | Gonapodya (strain JEL0774) |
| A042 | (lab id) NRRL 2543 | Rhizomucor pusillus NRRL 2543 | GCA_022702545.1 Phycomyces nitens | 0.9898 reads; 0.9986 assembly |  | ok (both) | READS+ASM | Phycomyces nitens |
| A043 | (lab id) RSA 1136 | Chaetocladium brefeldii RSA 1136 | GCA_017871115.1 Cunninghamella blakesleeana | 0.9864 reads; 0.9892 assembly |  | ok (both) | READS+ASM | Cunninghamella blakesleeana |

### Tier B

#### BFD genomes (57 left out of the reference from 2026.10.3 and 15 more from 2026.10.5; the others kept, flag kept_in_reference or size_outlier) (15)

| id | accession | current organism (BFD/NCBI) | GmG best hit | ANI | size | flags | evidence | suggested |
|---|---|---|---|---|---|---|---|---|
| B001 | GCA_035926115.1 | Elasticomyces elasticus |  | 0.9998 |  | audit43;hub_n=5;cluster95_others=0 | HUB | Eurotiomycetes (e.g. Exophiala sideris) |
| B002 | GCA_017654785.1 | Phenoliferia glacialis | GCA_017654825.1 Cryptococcus sp. T15_22C | 0.9951 |  | audit43;cluster95_others=1 | LOO+CLUSTER | Tremellomycetes (e.g. Cryptococcus sp. T15_22C) |
| B003 | GCA_017654825.1 | Cryptococcus sp. T15_22C | GCA_017654785.1 Phenoliferia glacialis | 0.9951 |  | audit43;cluster95_others=1 | LOO+CLUSTER | Microbotryomycetes (e.g. Phenoliferia glacialis) |
| B004 | GCA_002153405.1 | Mycosphaerella populi | GCA_039702135.1 Penicillium sp. 32 TS-2023 | 0.9930 |  | audit43;cluster95_others=1 | LOO+CLUSTER | Eurotiomycetes (e.g. Penicillium sp. 32 TS-2023) |
| B005 | GCA_039702135.1 | Penicillium sp. 32 TS-2023 | GCA_002153405.1 Mycosphaerella populi | 0.9930 |  | audit43;cluster95_others=1 | LOO+CLUSTER | Dothideomycetes (e.g. Mycosphaerella populi) |
| B006 | GCA_002994575.1 | Uromyces transversalis | GCA_002072375.1 Penicillium steckii | 0.9588 |  | audit43;cluster95_others=0;neighbour_in_reference | LOO | Penicillium steckii |
| B007 | GCA_022813125.1 | Septoriella neoarundinis | GCA_022814475.1 Cladosporium endophyticum | 0.9999 |  | kept_in_reference;cluster95_others=2;neighbour_in_reference | LOO+CLUSTER | Cladosporiales (e.g. Cladosporium endophyticum) |
| B008 | GCA_022814475.1 | Cladosporium endophyticum | GCA_022813125.1 Septoriella neoarundinis | 0.9999 |  | kept_in_reference;cluster95_others=2;neighbour_in_reference | LOO+CLUSTER | Pleosporales (e.g. Septoriella neoarundinis) |
| B009 | GCA_035927975.1 | Rachicladosporium monterosium | GCA_030553495.1 Friedmanniomyces endolithicus | 0.9994 |  | kept_in_reference;hub_n=12;cluster95_others=1;neighbour_in_reference | LOO+HUB+CLUSTER | Mycosphaerellales (e.g. Friedmanniomyces endolithicus) |
| B010 | GCA_030553495.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9994 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Cladosporiales (e.g. Rachicladosporium monterosium) |
| B011 | GCA_028028555.1 | Pseudopestalotiopsis chinensis | GCF_000516985.1 Pestalotiopsis fici W106-1 | 0.9913 |  | truthset14;cluster95_others=2 | REFSEQ+CLUSTER | Pestalotiopsis |
| B012 | GCA_051527695.1 | Neonectria sp. MB3469 | GCF_001636725.1 Cordyceps fumosorosea ARSEF 2679 | 0.9913 |  | truthset14;cluster95_others=1 | REFSEQ+CLUSTER | Cordyceps |
| B013 | GCA_003314395.1 | Albatrellus ellisii | GCA_003314755.1 Pleurotus eryngii | 0.9883 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Agaricales (e.g. Pleurotus eryngii) |
| B014 | GCA_003314755.1 | Pleurotus eryngii | GCA_003314395.1 Albatrellus ellisii | 0.9883 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Russulales (e.g. Albatrellus ellisii) |
| B015 | GCA_056321185.1 | Aspergillus destruens | GCF_028827035.1 Penicillium chrysogenum | 0.9836 |  | truthset14;cluster95_others=1 | REFSEQ+CLUSTER | Penicillium |

#### Lab or provider labels of read sets and assemblies, the SRA run, the decoy (10)

| id | accession | current organism (BFD/NCBI) | GmG best hit | ANI | size | flags | evidence | suggested |
|---|---|---|---|---|---|---|---|---|
| B016 | ERR14883901 / SAMEA118140668 / GCA_982397225.1 | Coccidioides immitis | GCA_982397405.1 Ophidiomyces ophidiicola (T55-F) | 0.9473 (reads, own genome hidden); 0.9503 (assembly, leave-one-out) |  | 5 other C. immitis genomes in the reference; next hit 19 shared hashes vs 772 | SRA+LOO | not C. immitis; Onygenaceae, nearest Ophidiomyces ophidiicola (ANI 0.95 is below species level, so the genus is not certain) |
| B017 | (lab id) NRRL 13127 | Umbelopsis ovata NRRL 13127 | GCA_054127015.1 Rhizopus delemar (Mucoromycetes; label is Umbelopsidomycetes) | 0.9976 |  | ok | ASM | Rhizopus delemar |
| B018 | (lab id) STR0027783 | Magnaporthe sp. QU1123 (IndigoAg) | GCA_001599375.1 Phoma herbarum | 0.9921 |  | ok | ASM | Phoma herbarum |
| B019 | (lab id) STR0027683 | Trichoderma sp. TK1246 (IndigoAg, provider 2017 name) | GCA_060370295.1 Aspergillus fumigatus | 0.9893 |  | ok | ASM | Aspergillus fumigatus |
| B020 | (lab id) STR0027657 | Magnaporthe sp. P3064 (IndigoAg) | GCA_027942835.1 Colletotrichum tropicale | 0.9884 |  | ok | ASM | Colletotrichum tropicale |
| B021 | (lab id) STR0027726 | Alternaria sp. DC477 (IndigoAg) | GCA_023625515.1 Aspergillus terreus | 0.9864 |  | ok | ASM | Aspergillus terreus |
| B022 | (lab id) STR0027671 | Lophiostoma (IndigoAg) | GCA_028829675.1 Penicillium tannophilum | 0.9864 |  | ok | ASM | Penicillium tannophilum |
| B023 | (lab id) JEL0847 | Paraphysoderma sedebokerense JEL0847 | GCA_027604345.1 Boothiomyces sp. JEL0866 (Chytridiomycota) | 0.9840 |  | ok; 1332 shared hashes | READS | Boothiomyces |
| B024 | (lab id) NRRL 5842 | Linnemannia nantahalensis NRRL 5842 | GCA_027405735.1 Mucor racemosus (Mucoromycota; label is Mortierellomycota) | 0.9631 |  | ok | ASM | Mucor |
| B025 | (lab id) NRRL 3022 | Conidiobolus incongruus NRRL 3022 | GCA_054601915.1 Gilbertella persicaria (Mucoromycota) | 0.9332 reads; 0.9391 assembly |  | ok (both) | READS+ASM | Gilbertella (Mucorales) or contaminated |

### Tier C

#### BFD genomes (57 left out of the reference from 2026.10.3 and 15 more from 2026.10.5; the others kept, flag kept_in_reference or size_outlier) (158)

| id | accession | current organism (BFD/NCBI) | GmG best hit | ANI | size | flags | evidence | suggested |
|---|---|---|---|---|---|---|---|---|
| C001 | GCA_024606205.1 | Alternaria sp. SPS-2 | GCA_003069565.1 Aspergillus fumigatus | 1.0000 |  | audit43;cluster95_others=209 | LOO+CLUSTER | Eurotiomycetes (e.g. Aspergillus fumigatus) |
| C002 | GCA_029300795.1 | Saccharomycetales sp. | GCF_000209165.1 Scheffersomyces stipitis | 0.9999 |  | audit43;hub_n=4;cluster95_others=5;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Scheffersomyces |
| C003 | GCA_029290875.1 | Saccharomycetales sp. | GCA_000755205.1 Meyerozyma caribbica | 0.9987 |  | audit43;hub_n=18;cluster95_others=18;neighbour_in_reference | LOO+HUB+CLUSTER | Meyerozyma guilliermondii (17 of 18) |
| C004 | GCA_029301555.1 | Saccharomycetales sp. | GCA_025531575.1 Pichia sp. PH4023 | 0.9982 |  | audit43;hub_n=6;cluster95_others=3;neighbour_in_reference | LOO+HUB+REFSEQ+CLUSTER | Pichia (3 of 3; e.g. Pichia sp. PH4023) |
| C005 | GCA_029301665.1 | Saccharomycetales sp. | GCA_025531575.1 Pichia sp. PH4023 | 0.9971 |  | audit43;cluster95_others=0;neighbour_in_reference | LOO | Pichia sp. PH4023 |
| C006 | GCA_029301595.1 | Saccharomycetales sp. | GCA_012295375.1 Brettanomyces bruxellensis | 0.9958 |  | audit43;cluster95_others=0;neighbour_in_reference | LOO | Brettanomyces bruxellensis |
| C007 | GCA_029295605.1 | Saccharomycetales sp. | GCA_023623335.1 Meyerozyma caribbica | 0.9940 |  | audit43;hub_n=9;cluster95_others=9;neighbour_in_reference | LOO+HUB+CLUSTER | Meyerozyma caribbica (5 of 9) |
| C008 | GCA_057557315.1 | Saccharomycetes sp. | GCF_011074885.1 Brettanomyces bruxellensis | 0.9922 |  | audit43;hub_n=10;cluster95_others=10 | LOO+HUB+REFSEQ+CLUSTER | Brettanomyces bruxellensis (9 of 10) |
| C009 | GCF_011074885.1 | Brettanomyces bruxellensis | GCA_057557315.1 Saccharomycetes sp. | 0.9922 |  | audit43;cluster95_others=10 | LOO+CLUSTER | Saccharomycetes (e.g. Saccharomycetes sp.) |
| C010 | GCA_037042185.1 | Sordariomycetes sp. S80-01Y | GCA_037074925.2 Tetracladium maxilliforme | 0.9900 |  | audit43;cluster95_others=1 | LOO+CLUSTER | Leotiomycetes (e.g. Tetracladium maxilliforme) |
| C011 | GCA_037074925.2 | Tetracladium maxilliforme | GCA_037042185.1 Sordariomycetes sp. S80-01Y | 0.9900 |  | audit43;cluster95_others=1 | LOO+CLUSTER | Sordariomycetes (e.g. Sordariomycetes sp. S80-01Y) |
| C012 | GCA_016906875.1 | Brettanomyces sp. HC-2020a | GCA_057557315.1 Saccharomycetes sp. | 0.9861 |  | audit43;cluster95_others=10 | LOO+CLUSTER | Saccharomycetes sp. |
| C013 | GCA_052696405.1 | Microbotryomycetes sp. NB124-2 | GCA_054471015.1 Saccharomycotina sp. YJ-2025b | 0.9819 |  | audit43;cluster95_others=1 | LOO+CLUSTER | Ascomycota (e.g. Saccharomycotina sp. YJ-2025b) |
| C014 | GCA_054471015.1 | Saccharomycotina sp. YJ-2025b | GCA_052696405.1 Microbotryomycetes sp. NB124-2 | 0.9819 |  | audit43;cluster95_others=1 | LOO+CLUSTER | Basidiomycota (e.g. Microbotryomycetes sp. NB124-2) |
| C015 | GCA_057978885.1 | Chaetomiaceae sp. 537B | GCA_056320565.1 Penicillium chrysogenum | 0.9776 |  | audit43;cluster95_others=118;neighbour_in_reference | LOO+CLUSTER | Penicillium chrysogenum (89 of 118) |
| C016 | GCA_001599655.1 | Golubevia pallescens | GCA_030280995.1 Hemileia vastatrix | 0.9699 |  | audit43;cluster95_others=0;neighbour_in_reference | LOO | Hemileia vastatrix |
| C017 | GCA_057383725.1 | Candida tropicalis |  | 0.9545 |  | audit43;hub_n=19;cluster95_others=111 | HUB+CLUSTER | Saccharomycetes (e.g. Nakaseomyces glabratus) |
| C018 | GCA_040938115.1 | Rigidonotus glomeratus | GCA_037892935.1 Phlebopus portentosus | 1.0000 |  | kept_in_reference;cluster95_others=5;neighbour_in_reference | LOO+CLUSTER | Boletales (e.g. Phlebopus portentosus) |
| C019 | GCA_023087545.1 | Lasiodiplodia pseudotheobromae | GCA_023087585.1 Bipolaris maydis | 0.9995 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pleosporales (e.g. Bipolaris maydis) |
| C020 | GCA_937900885.1 | uncultured Penicillium | GCF_028828025.1 Penicillium rubens | 0.9972 |  | truthset14;cluster95_others=118 | REFSEQ+CLUSTER | Penicillium |
| C021 | GCA_054785315.1 | Glaciozyma watsonii | GCA_937872655.1 Leucosporidium creatinivorum | 0.9960 |  | kept_in_reference;cluster95_others=3;neighbour_in_reference | LOO+CLUSTER | Leucosporidiales (e.g. Leucosporidium creatinivorum) |
| C022 | GCA_054785265.1 | Glaciozyma watsonii | GCA_937872655.1 Leucosporidium creatinivorum | 0.9950 |  | kept_in_reference;cluster95_others=3;neighbour_in_reference | LOO+CLUSTER | Leucosporidium creatinivorum |
| C023 | GCA_054785335.1 | Glaciozyma watsonii | GCA_937872655.1 Leucosporidium creatinivorum | 0.9949 |  | kept_in_reference;cluster95_others=3;neighbour_in_reference | LOO+CLUSTER | Leucosporidium creatinivorum |
| C024 | GCA_018153875.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9938 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Mycosphaerellales (e.g. Pseudocercospora fijiensis) |
| C025 | GCA_029169105.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9926 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C026 | GCA_002111425.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9923 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C027 | GCA_029169125.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9922 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C028 | GCA_023086975.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9922 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C029 | GCA_029169085.1 | Lasiodiplodia pseudotheobromae | GCA_023079265.1 Alternaria brassicicola | 0.9921 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pleosporales (e.g. Alternaria brassicicola) |
| C030 | GCA_023087325.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9917 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C031 | GCA_060251805.1 | Wickerhamomyces sp. 50C-3 | GCA_030574675.1 Candida sp. | 0.9908 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Saccharomycetales (e.g. Candida sp.) |
| C032 | GCA_030574675.1 | Candida sp. | GCA_060251805.1 Wickerhamomyces sp. 50C-3 | 0.9908 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Phaffomycetales (e.g. Wickerhamomyces sp. 50C-3) |
| C033 | GCA_022814635.1 | Cladosporium colombiae | GCA_022813475.1 Pseudopyrenochaeta terrestris | 0.9907 |  | kept_in_reference;cluster95_others=3;neighbour_in_reference | LOO+CLUSTER | Pleosporales (e.g. Pseudopyrenochaeta terrestris) |
| C034 | GCA_030560805.1 | Suzukiozyma palmioleophila | GCA_937898345.1 uncultured Saccharomycetales | 0.9901 |  | kept_in_reference;cluster95_others=4;neighbour_in_reference | LOO+CLUSTER | uncultured Saccharomycetales |
| C035 | GCA_049997415.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9900 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C036 | GCA_029295655.1 | Hypocreales sp. | GCF_000392275.1 Phaeoacremonium minimum | 0.9895 |  | kept_in_reference;cluster95_others=0;neighbour_in_reference | LOO | Phaeoacremonium minimum |
| C037 | GCA_009829825.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9895 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C038 | GCA_977013645.2 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9892 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C039 | GCF_012971845.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9888 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C040 | GCA_023087405.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9887 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C041 | GCA_025602875.1 | Suzukiozyma palmioleophila | GCA_937898345.1 uncultured Saccharomycetales | 0.9886 |  | kept_in_reference;cluster95_others=4;neighbour_in_reference | LOO+CLUSTER | uncultured Saccharomycetales |
| C042 | GCA_977013665.2 | Lasiodiplodia brasiliensis | GCA_011057605.1 Pseudocercospora fijiensis | 0.9886 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C043 | GCA_008931885.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9883 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C044 | GCA_002901145.1 | Cladosporium cladosporioides | GCA_022813475.1 Pseudopyrenochaeta terrestris | 0.9835 |  | kept_in_reference;cluster95_others=3;neighbour_in_reference | LOO+CLUSTER | Pseudopyrenochaeta terrestris |
| C045 | GCA_024195665.1 | Pleurotus ostreatus | GCA_017499575.1 Sanguinoderma rugosum | 0.9820 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Polyporales (e.g. Sanguinoderma rugosum) |
| C046 | GCA_019775475.1 | Cladosporium colombiae | GCA_022813475.1 Pseudopyrenochaeta terrestris | 0.9803 |  | kept_in_reference;cluster95_others=3;neighbour_in_reference | LOO+CLUSTER | Pseudopyrenochaeta terrestris |
| C047 | GCA_023373445.1 | Lasiodiplodia theobromae | GCA_023079265.1 Alternaria brassicicola | 0.9803 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Alternaria brassicicola |
| C048 | GCA_001956935.2 | Pleurotus ostreatus | GCA_017499575.1 Sanguinoderma rugosum | 0.9798 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Sanguinoderma rugosum |
| C049 | GCA_947034855.1 | Pleurotus ostreatus | GCA_017499575.1 Sanguinoderma rugosum | 0.9781 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Sanguinoderma rugosum |
| C050 | GCA_056149315.1 | Pleurotus ostreatus | GCA_017499575.1 Sanguinoderma rugosum | 0.9748 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Sanguinoderma rugosum |
| C051 | GCA_056149245.1 | Pleurotus ostreatus | GCA_017499575.1 Sanguinoderma rugosum | 0.9745 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Sanguinoderma rugosum |
| C052 | GCA_051307915.1 | Pleurotus ostreatus | GCA_017499575.1 Sanguinoderma rugosum | 0.9745 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Sanguinoderma rugosum |
| C053 | GCA_030512325.1 | Lasiodiplodia hormozganensis | GCA_011057605.1 Pseudocercospora fijiensis | 0.9713 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C054 | GCA_008079295.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9713 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C055 | GCA_019677325.2 | Pleurotus cornucopiae | GCA_017499575.1 Sanguinoderma rugosum | 0.9711 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Sanguinoderma rugosum |
| C056 | GCA_042257065.1 | Naganishia uzbekistanensis | GCA_900079645.1 Cryptococcus sp. 05_00 | 0.9702 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Tremellales (e.g. Cryptococcus sp. 05_00) |
| C057 | GCA_900079645.1 | Cryptococcus sp. 05_00 | GCA_042257065.1 Naganishia uzbekistanensis | 0.9702 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Filobasidiales (e.g. Naganishia uzbekistanensis) |
| C058 | GCA_019395345.1 | Pleurotus floridanus | GCA_017499575.1 Sanguinoderma rugosum | 0.9700 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Sanguinoderma rugosum |
| C059 | GCA_003313735.1 | Pleurotus platypus | GCA_017499575.1 Sanguinoderma rugosum | 0.9689 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Sanguinoderma rugosum |
| C060 | GCA_002778095.1 | Fragosphaeria purpurea | GCA_001636815.1 Niveomyces insectorum | 0.9687 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Hypocreales (e.g. Niveomyces insectorum) |
| C061 | GCA_001636815.1 | Niveomyces insectorum | GCA_002778095.1 Fragosphaeria purpurea | 0.9687 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Ophiostomatales (e.g. Fragosphaeria purpurea) |
| C062 | GCA_022264855.1 | Lenzites betulinus | GCA_059750515.1 Cortinarius caperatus | 0.9678 |  | kept_in_reference;cluster95_others=0;neighbour_in_reference | LOO | Cortinarius caperatus |
| C063 | GCA_030553335.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9634 |  | kept_in_reference;cluster95_others=11;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C064 | GCA_046563505.1 | Gemmamyces piceae | GCA_046127395.1 Dichomera gemmicola | 0.9620 |  | kept_in_reference;cluster95_others=6;neighbour_in_reference | LOO+CLUSTER | Dichomera gemmicola |
| C065 | GCA_046563545.1 | Gemmamyces piceae | GCA_046127395.1 Dichomera gemmicola | 0.9618 |  | kept_in_reference;cluster95_others=6;neighbour_in_reference | LOO+CLUSTER | Dichomera gemmicola |
| C066 | GCA_041146655.1 | Gemmamyces piceae | GCA_046127395.1 Dichomera gemmicola | 0.9616 |  | kept_in_reference;cluster95_others=6;neighbour_in_reference | LOO+CLUSTER | Dichomera gemmicola |
| C067 | GCA_003314385.1 | Gomphus sp. MG54 | GCA_003314275.1 Chroogomphus rutilus | 0.9606 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Boletales (e.g. Chroogomphus rutilus) |
| C068 | GCA_003314275.1 | Chroogomphus rutilus | GCA_003314385.1 Gomphus sp. MG54 | 0.9606 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Gomphales (e.g. Gomphus sp. MG54) |
| C069 | GCA_041430495.1 | Kazachstania sp. VP_05 | GCF_036370985.1 Arxiozyma heterogenica | 0.9605 |  | truthset14;cluster95_others=13 | REFSEQ+CLUSTER | Arxiozyma |
| C070 | GCA_027569505.1 | Chloridium humicola | GCA_004802645.1 Amesomyces nigricolor | 0.9602 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Sordariales (e.g. Amesomyces nigricolor) |
| C071 | GCA_004802645.1 | Amesomyces nigricolor | GCA_027569505.1 Chloridium humicola | 0.9602 |  | kept_in_reference;cluster95_others=1;neighbour_in_reference | LOO+CLUSTER | Chaetosphaeriales (e.g. Chloridium humicola) |
| C072 | GCA_023089105.1 | Lasiodiplodia citricola | GCA_011057605.1 Pseudocercospora fijiensis | 0.9602 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C073 | GCA_035621335.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9597 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C074 | GCA_041146615.1 | Gemmamyces piceae | GCA_046127395.1 Dichomera gemmicola | 0.9594 |  | kept_in_reference;cluster95_others=6;neighbour_in_reference | LOO+CLUSTER | Dichomera gemmicola |
| C075 | GCA_052857205.1 | Lasiodiplodia theobromae | GCA_011057605.1 Pseudocercospora fijiensis | 0.9589 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C076 | GCA_977013655.2 | Lasiodiplodia euphorbiaceicola | GCA_011057605.1 Pseudocercospora fijiensis | 0.9586 |  | kept_in_reference;cluster95_others=43;neighbour_in_reference | LOO+CLUSTER | Pseudocercospora fijiensis |
| C077 | GCA_009666835.1 | Saccharomyces cerevisiae x Saccharomyces eubayanus x Saccharomyces uvarum | GCA_000350065.2 Colletotrichum orbiculare | 0.9573 |  | kept_in_reference;cluster95_others=0;neighbour_in_reference | LOO | Colletotrichum orbiculare |
| C078 | GCA_030926145.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9550 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C079 | GCA_016618145.1 | Sanghuangporus lonicericola | GCA_001792895.1 Lentinus brumalis | 0.9534 |  | kept_in_reference;cluster95_others=0;neighbour_in_reference | LOO | Lentinus brumalis |
| C080 | GCA_030553205.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9519 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C081 | GCA_030553215.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9518 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C082 | GCA_030553295.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9515 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C083 | GCA_030553275.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9514 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C084 | GCA_030553195.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9514 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C085 | GCA_041146935.1 | Gemmamyces piceae | GCA_046127395.1 Dichomera gemmicola | 0.9514 |  | kept_in_reference;cluster95_others=6;neighbour_in_reference | LOO+CLUSTER | Dichomera gemmicola |
| C086 | GCA_030553175.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9514 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C087 | GCA_030412055.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9514 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C088 | GCA_030553285.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9513 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C089 | GCA_030553475.1 | Friedmanniomyces endolithicus | GCA_035927975.1 Rachicladosporium monterosium | 0.9511 |  | kept_in_reference;cluster95_others=9;neighbour_in_reference | LOO+CLUSTER | Rachicladosporium monterosium |
| C090 | GCA_000167075.2 | Saccharomyces kudriavzevii |  |  | 1.75 | size_outlier | SIZE |  |
| C091 | GCA_000257045.1 | Saccharomyces kudriavzevii |  |  | 1.6 | size_outlier | SIZE |  |
| C092 | GCA_000697275.1 | Rhizopus microsporus |  |  | 2.73 | size_outlier | SIZE |  |
| C093 | GCA_000961905.2 | Agroathelia rolfsii |  |  | 1.69 | size_outlier | SIZE |  |
| C094 | GCA_001950435.1 | Zygosaccharomyces rouxii |  |  | 1.7 | size_outlier | SIZE |  |
| C095 | GCA_001972345.1 | Zygosaccharomyces rouxii |  |  | 1.69 | size_outlier | SIZE |  |
| C096 | GCA_002104975.1 | Neocallimastix californiae |  |  | 1.68 | size_outlier | SIZE |  |
| C097 | GCA_002196765.1 | Millerozyma farinosa |  |  | 1.9 | size_outlier | SIZE |  |
| C098 | GCA_002317995.1 | Diplocarpon rosae |  |  | 1.72 | size_outlier | SIZE |  |
| C099 | GCA_003314755.1 | Pleurotus eryngii |  |  | 2.33 | size_outlier | SIZE |  |
| C100 | GCA_003327635.1 | Saccharomyces kudriavzevii |  |  | 2.01 | size_outlier | SIZE |  |
| C101 | GCA_013181225.1 | Saccharomyces cerevisiae |  |  | 1.76 | size_outlier | SIZE |  |
| C102 | GCA_013181245.1 | Saccharomyces cerevisiae |  |  | 1.64 | size_outlier | SIZE |  |
| C103 | GCA_013183865.1 | Saccharomyces cerevisiae |  |  | 1.89 | size_outlier | SIZE |  |
| C104 | GCA_013183875.1 | Saccharomyces cerevisiae |  |  | 1.86 | size_outlier | SIZE |  |
| C105 | GCA_016808315.1 | Rhodotorula toruloides |  |  | 1.81 | size_outlier | SIZE |  |
| C106 | GCA_018343895.1 | Agroathelia rolfsii |  |  | 1.7 | size_outlier | SIZE |  |
| C107 | GCA_019827475.1 | Penicillium cf. chrysogenum |  |  | 2.22 | size_outlier | SIZE |  |
| C108 | GCA_020521415.1 | Zygosaccharomyces rouxii |  |  | 1.73 | size_outlier | SIZE |  |
| C109 | GCA_024195665.1 | Pleurotus ostreatus |  |  | 1.75 | size_outlier | SIZE |  |
| C110 | GCA_024256405.1 | Debaryomyces hansenii |  |  | 1.8 | size_outlier | SIZE |  |
| C111 | GCA_024761585.1 | Hortaea werneckii |  |  | 1.61 | size_outlier | SIZE |  |
| C112 | GCA_029035315.1 | Colletotrichum liriopes |  |  | 1.66 | size_outlier | SIZE |  |
| C113 | GCA_029292925.1 | Ophiocordycipitaceae sp. |  |  | 3.31 | size_outlier | SIZE |  |
| C114 | GCA_029292965.1 | Ophiocordycipitaceae sp. |  |  | 2.36 | size_outlier | SIZE |  |
| C115 | GCA_029292995.1 | Ophiocordycipitaceae sp. |  |  | 3.23 | size_outlier | SIZE |  |
| C116 | GCA_029294875.1 | Ascomycota sp. |  |  | 2.02 | size_outlier | SIZE |  |
| C117 | GCA_029295915.1 | Ascomycota sp. |  |  | 2.25 | size_outlier | SIZE |  |
| C118 | GCA_029300675.1 | Saccharomycetaceae sp. |  |  | 3.21 | size_outlier | SIZE |  |
| C119 | GCA_029303725.1 | Saccharomycetaceae sp. |  |  | 3.15 | size_outlier | SIZE |  |
| C120 | GCA_029303885.1 | Saccharomycetaceae sp. |  |  | 1.89 | size_outlier | SIZE |  |
| C121 | GCA_030913545.1 | Cryptococcus neoformans |  |  | 1.8 | size_outlier | SIZE |  |
| C122 | GCA_041430195.1 | Agaricus bisporus |  |  | 2.18 | size_outlier | SIZE |  |
| C123 | GCA_043234885.1 | Exophiala dermatitidis |  |  | 1.72 | size_outlier | SIZE |  |
| C124 | GCA_043234895.1 | Exophiala dermatitidis |  |  | 1.72 | size_outlier | SIZE |  |
| C125 | GCA_054960005.1 | Buchwaldoboletus xylophilus |  |  | 2.65 | size_outlier | SIZE |  |
| C126 | GCA_055773995.1 | Malassezia globosa |  |  | 2.12 | size_outlier | SIZE |  |
| C127 | GCA_055774135.1 | Sphaerulina musiva |  |  | 1.89 | size_outlier | SIZE |  |
| C128 | GCA_058127325.1 | Buchwaldoboletus xylophilus |  |  | 2.79 | size_outlier | SIZE |  |
| C129 | GCA_058127345.1 | Buchwaldoboletus xylophilus |  |  | 2.8 | size_outlier | SIZE |  |
| C130 | GCA_058127365.1 | Buchwaldoboletus xylophilus |  |  | 2.77 | size_outlier | SIZE |  |
| C131 | GCA_058127385.1 | Buchwaldoboletus xylophilus |  |  | 2.62 | size_outlier | SIZE |  |
| C132 | GCA_058127405.1 | Buchwaldoboletus xylophilus |  |  | 2.75 | size_outlier | SIZE |  |
| C133 | GCA_058127425.1 | Buchwaldoboletus xylophilus |  |  | 2.86 | size_outlier | SIZE |  |
| C134 | GCA_058129735.1 | Buchwaldoboletus xylophilus |  |  | 2.8 | size_outlier | SIZE |  |
| C135 | GCA_900411385.1 | Zygosaccharomyces rouxii |  |  | 1.81 | size_outlier | SIZE |  |
| C136 | GCA_900682665.1 | Saccharomyces kudriavzevii |  |  | 1.94 | size_outlier | SIZE |  |
| C137 | GCA_905219595.1 | Rhizoctonia solani |  |  | 1.68 | size_outlier | SIZE |  |
| C138 | GCA_918252685.1 | Saccharomyces kudriavzevii |  |  | 1.99 | size_outlier | SIZE |  |
| C139 | GCA_918282095.1 | Saccharomyces kudriavzevii |  |  | 1.98 | size_outlier | SIZE |  |
| C140 | GCA_979262465.1 | Zygosaccharomyces siamensis |  |  | 1.95 | size_outlier | SIZE |  |
| C141 | GCA_979262545.1 | Zygosaccharomyces rouxii |  |  | 1.69 | size_outlier | SIZE |  |
| C142 | GCA_979262605.1 | Torulaspora microellipsoides |  |  | 1.9 | size_outlier | SIZE |  |
| C143 | GCA_979262865.1 | Zygosaccharomyces rouxii |  |  | 1.7 | size_outlier | SIZE |  |
| C144 | GCA_979263555.1 | Zygosaccharomyces siamensis |  |  | 1.95 | size_outlier | SIZE |  |
| C145 | GCA_979719475.1 | Torulaspora microellipsoides |  |  | 1.87 | size_outlier | SIZE |  |
| C146 | GCA_980284945.1 | Zygosaccharomyces rouxii |  |  | 1.69 | size_outlier | SIZE |  |
| C147 | GCA_980285225.1 | Zygosaccharomyces rouxii |  |  | 1.61 | size_outlier | SIZE |  |
| C148 | GCA_980286765.1 | Zygosaccharomyces rouxii |  |  | 1.63 | size_outlier | SIZE |  |
| C149 | GCA_980286925.1 | Zygosaccharomyces rouxii |  |  | 1.64 | size_outlier | SIZE |  |
| C150 | GCA_980287995.1 | Zygosaccharomyces siamensis |  |  | 1.95 | size_outlier | SIZE |  |
| C151 | GCA_980288055.1 | Zygosaccharomyces siamensis |  |  | 1.95 | size_outlier | SIZE |  |
| C152 | GCA_980288285.1 | Zygosaccharomyces rouxii |  |  | 1.69 | size_outlier | SIZE |  |
| C153 | GCA_980290965.1 | Zygosaccharomyces rouxii |  |  | 1.68 | size_outlier | SIZE |  |
| C154 | GCA_980291795.1 | Zygosaccharomyces rouxii |  |  | 1.68 | size_outlier | SIZE |  |
| C155 | GCA_980292995.1 | Zygosaccharomyces rouxii |  |  | 1.65 | size_outlier | SIZE |  |
| C156 | GCA_980294845.1 | Zygosaccharomyces rouxii |  |  | 1.68 | size_outlier | SIZE |  |
| C157 | GCA_980695105.1 | Zygosaccharomyces rouxii |  |  | 1.67 | size_outlier | SIZE |  |
| C158 | GCF_947243775.1 | Saccharomyces kudriavzevii |  |  | 2.0 | size_outlier | SIZE |  |

#### Lab or provider labels of read sets and assemblies, the SRA run, the decoy (10)

| id | accession | current organism (BFD/NCBI) | GmG best hit | ANI | size | flags | evidence | suggested |
|---|---|---|---|---|---|---|---|---|
| C159 | (lab id) Cryphonectria__CPWV1912 | Cryphonectria parasitica CPWV1912 (KassonLab) | n/a Neonectria punicea | 0.9972 |  | ok; input is archived herbarium DNA, a metagenome | ASM |  |
| C160 | (lab id) NRRL 22899 | Circinomucor circinelloides NRRL 22899 (reads sheet) / Mucor circinelloides (assembly sheet) | GCA_000401635.1 Mucor circinelloides | 0.9969 reads; 0.9992 assembly |  | ok (both) | READS+ASM | Mucor circinelloides (or Circinomucor circinelloides) |
| C161 | (lab id) STR0027688 | Neosetophoma samarorum (IndigoAg) | GCA_023628855.1 Coniochaeta sp. PK30-01 | 0.9781 |  | ok | ASM | Coniochaeta |
| C162 | (lab id) CBS 441.59 | Syncephalastrum racemosum CBS 441.59 | GCA_017871115.1 Cunninghamella blakesleeana (reads); Syncephalastrum racemosum (assembly) | 0.9736 reads; 0.9654 assembly |  | ok (both) | READS vs ASM |  |
| C163 | (lab id) STR0027763 | Neurospora sp. PS0169 (IndigoAg) | GCA_023626415.1 Penicillium sp. S1701 | 0.9714 |  | ok | ASM | Penicillium |
| C164 | (lab id) STR0027680 | Diaporthe sp. P3581 (IndigoAg) | GCA_022385635.1 Xylaria scruposa | 0.9623 |  | ok | ASM | Xylaria |
| C165 | (lab id) OphiocordycepsSymbiont__HPPF3 | Ophiocordyceps sp. HPPF3 (KassonLab) | n/a Massospora platypediae (Entomophthoromycota) | 0.9580 |  | ok; fragmented; input is a metagenome | ASM |  |
| C166 | (lab id) STR0027681 | Nemania sp. QU1039 (IndigoAg) | GCA_001879275.1 Coniochaeta ligniaria | 0.9548 |  | ok | ASM | Coniochaeta |
| C167 | (lab id) NRRL A-11121 | Mucor mucedo NRRL A-11121 | GCA_001986395.2 Irpex lacteus (Basidiomycota) | 0.9390 |  | ok; 161 shared of 604 hashes (small sketch) | READS | unknown (Irpex-like contaminant) |
| C168 | GCF_965202585.1 | Thalia democratica (salp) | 9 ZyGoLife assemblies + 1 read set Conidiobolus/Neoconidiobolus (7), Mucor hiemalis (2), Conidiobolus antarcticus ARSEF 691 reads | 0.83-0.913 (best hit of 9 assemblies, 1 read set) |  | false not_fungal on Entomophthoromycota (2 assemblies at ANI 0.9001/0.9055, 1 read set at 0.913) | ASM+READS | unknown; possible fungal contaminant contigs in the assembly, or a conserved repeat |

<!-- /TABLES -->

## Ambiguities, and what was not decided

* **Exclusion status** is taken from `reports/label_audit/reference_exclusions.txt` (57 ids, all found in `data/external/bfd_meta.tsv`); the packed index files were not opened. The published `index-2026.10.3` was built before the dangling-cluster fix (`reports/validation/clusters_vs_skani.md`) but with the same 57 exclusions.
* **The 14 truth-set exclusions** have no per-genome evidence file in the repository; only 4 have a hand review (`truth_set_reviews.tsv`). Their rows are rebuilt from the first truth-set run and from the clusters, so they may differ in detail from what #90 had in mind when it chose the 14.
* **Which side of a pair is wrong** is not known for the B pairs (e.g. *Phenoliferia glacialis* GCA_017654785.1 against *Cryptococcus* sp. GCA_017654825.1, both left out). For a pair where one side was left out and the other kept, the kept side is listed with its own row if it is itself odd.
* **Pseudopestalotiopsis**: GCA_028028555.1 (left out) is the probable error against the RefSeq *Pestalotiopsis fici*, but the *Pseudopestalotiopsis theae* genome GCA_015881745.1 that stays in the reference sits at ANI 0.990 from the same RefSeq genome, so a genus-level naming question remains for the whole group (not decided).
* **Order-level taxonomy** (e.g. Mycosphaerellales/Cladosporiales/Capnodiales, Polyporales/Agaricales): 70 kept genomes have an order-level conflict that looks like a classification difference (tier C); some may be real mislabels, as the tier-A order-level cases show.
* **Not listed one by one**: 2,662 conflicts in `log_conflicts.tsv`, of which 2,490 (1,245 genomes) have one of the 57 as the conflicting neighbour and are resolved by the exclusion; the other 89 rows (85 genomes, all between kept genomes) are handled in the kept-genome rows. IndigoAg genus or family differences that look like renamed genera (14 at ANI >= 0.95, e.g. *Lecythophora* to *Coniochaeta*, *Kretzschmaria* to *Ustulina*) and the many ZyGoLife genus differences were not reviewed one by one: only strains with reads and assembly, and class-level changes, are listed.
* **External sets without names**: the 266 Antarctic MAGs and the CCFEE isolates have no names to dispute (`external_antarctica.md`, `external_ccfee.md`), so none is listed; the Kasson metagenome and herbarium inputs are listed only as tier C.
* **Provider names are not curated**: the IndigoAg names are the provider's 2017 strain names (genus guessed); the top hit can itself be a mislabelled reference genome (it was not one of the 57 for any lab row, checked) so these B rows may be the reference's fault.

## Companion file

`reports/label_audit/upstream_report.tsv`: the same rows, one per genome or read set, columns `tier, row_id, kind, accession, bfd_id, current_organism, current_lineage, best_hit_accession, best_hit_taxon, ani, size_ratio, flags, evidence_type, evidence, suggested_taxon, suggested_action, in_reference_2026.10.3, sources`.
