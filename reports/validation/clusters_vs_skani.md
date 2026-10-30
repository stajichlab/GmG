# The 95% ANI clusters against skani (#47)

The `cluster` column of the index (greedy clustering of 23,669 BFD genomes into 6,104 clusters at containment ANI >= 0.95, k = 21, scaled 10000, centroids in
the BFD priority order) was checked against **skani 0.3.2** (default settings, `triangle -E`; skani reports ANI from alignments and is in its "learned ANI"
mode, which adjusts the value with a model trained on MAGs). Genomes are the BFD clean assemblies. Per-pair output stays on the HPCC
(`GmG_runs/cluster_skani`); `scripts/cluster_vs_skani.py` and `scripts/hpcc/27_cluster_skani.sbatch` redo it (a different `--seed` gives another sample).
Sample: 124 random clusters of size >= 2 (25 per size class 2, 3-5, 6-20, 21-100, > 100; up to 30 members), the four largest clusters of each of eight
species-complex genera (up to 60 members, all pairs), and all centroids of those genera (up to 150).

## Findings

1. **Members are mostly where gmg says, with a shift at the edge.** Of 1,661 random member-centroid pairs (gmg ANI >= 0.95 by construction) skani gives
   ANI < 95% for **91 (5.5%, 95% CI 4.5-6.7%)**, and nearly all of them sit at the boundary: **70 of the 85 pairs with gmg ANI 0.95-0.96 (82%)** have skani ANI below 0.95
   (mean 0.945), against 21 of 1,576 (1.3%) for gmg ANI >= 0.96 and **0 of 851 at >= 0.99**. gmg's containment ANI is on average 0.0028 higher than skani's
   (median 0.0015; 5th percentile -0.0115). So a "95% cluster" means about 94-95% by alignment at its edge, and clearly >= 95% inside.
   The *Fusarium fujikuroi* complex cluster (centroid *F. chuoi*, 84 members: *F. proliferatum*, *F. fujikuroi*, *F. annulatum*, ...) is the clearest case: members at gmg
   ANI 0.950-0.953 have skani ANI 0.932-0.942.
2. **Species complexes chain, as expected from a centroid rule.** Clusters that hold several named species of a complex have member pairs below 95%:
   *Cryptococcus gattii* complex (22 members, 6 species: 58% of pairs < 95%, 13 of 21 members < 95% of the centroid), a cluster of 83 *Aspergillus* with 12 species labels (35% of
   pairs), the *A. flavus/oryzae/parasiticus/sojae* cluster (518 members, 9.7%), *F. graminearum* clade (466, 6.6%), the *F. oxysporum* complex (810, 3.3%), *Saccharomyces*
   hybrid clusters (14-16%). Most large single-species clusters (*S. cerevisiae* 1,302, *C. auris* 950, *C. albicans* 134, *A. fumigatus* 416, *C. neoformans* 187) have 0-0.8% pairs below 95% and a
   minimum pair ANI of 0.947-0.991. Across the single-species clusters tested the minimum is 0.947-0.991, inside the complex clusters 0.925-0.94.
3. **Few clusters that skani would merge, and they are explained by assembly size.** Among all centroid pairs of the eight genera (23,822 pairs), skani gives
   >= 95% for 41 pairs (0.17%), 33 of them >= 97%: *Aspergillus* 5 of 8,515, *Fusarium* 19 of 7,260, *Saccharomyces* 14 of 136, *Candida* 2 of 136, *Penicillium* 1 of 7,626,
   none in *Candidozyma*, *Cryptococcus*, *Trichophyton*. They are same-species pairs (*S. kudriavzevii* / *S. kudriavzevii* 0.998, *C. albicans* / *C. albicans* 0.995,
   *P. chrysogenum* / *P. cf. chrysogenum* 0.998). Why gmg separates them: its ANI is the share of the *member's* k-mers found in the centroid; a member larger than the
   centroid (an incomplete centroid assembly, extra contigs, a different ploidy) scores low although the centroid is entirely inside the member. The sizes agree: the
   median size ratio (larger / smaller number of hashes) of the split pairs is **1.7-4.9**, of all pairs of the same genera 1.1-1.9.

## Consequences

* Species-level use of `cluster`: a cluster is a good proxy for a species where the species is well defined and a loose group in complexes (the report should
  not present a cluster as a species). This is the "species complex" caveat of `docs/limitations.md`, now with numbers.
* The edge shift is a calibration point, not a bug: if clusters are to mean "skani ANI >= 95%" the cutoff for the clustering would be about 0.955-0.96 in gmg ANI
  (at gmg ANI 0.96 only 4% of pairs are below 0.95 by skani). Not changed here: `cluster` is metadata, not used for the calls.
* Splits come from the asymmetry of the clustering containment. A cheap fix, if clusters ever need to be tighter: for clustering use the maximum of the two containments
  (sourmash "max containment"), or order centroids by completeness first. Not needed for calling.

## Cluster ids between releases

Clusters are computed once and carried unchanged: **2026.10, 2026.10.1 and 2026.10.2 have identical ids (23,669 of 23,669 genomes, adjusted Rand index 1.0000)**. Release 2026.10.3
removes the 57 ambiguous-label genomes (`reports/label_audit/reference_exclusions.txt`) and keeps the ids of the others: 23,612 of 23,612 genomes
keep their id, **but 15 of the 57 removed genomes were centroids, so 3 clusters (23 genomes) now carry the id of a genome that is not in the index** (8 more clusters vanish,
all their members were removed). `scripts/reseed_clusters.py` gives such a cluster a new centroid (the first remaining member in the BFD priority order);
`15_release_index.sbatch` now calls it in place of the plain exclusion filter. **The published `index-2026.10.3` was built before this and still has the 3 dangling ids**
(23 genomes; the `cluster` field of a result whose top hit is one of them names a genome that is not in the index; calls are not affected). The next release build fixes it; replacing the
2026.10.3 assets would change their checksums, which is a decision for the maintainers.

## Tables (generated by `scripts/cluster_vs_skani.py score`)

### 1. Members against their centroid (stratified random clusters)

124 clusters (sizes 2 to >100, up to 30 members each), 1661 member-centroid pairs. gmg's containment ANI (k=21, scaled 10000) is >= 0.95 for all of them by construction. **skani ANI < 95% (or not reported, i.e. < ~80% / too little alignment) for 91 of 1661 (5.5%; 95% CI 4.5% to 6.7%).**

Where skani reports a value (1661 pairs): mean skani - gmg = -0.0028, median -0.0015, 5th percentile -0.0115, 95th +0.0014. skani ANI < 0.95 for 91 of them.

| gmg ANI to centroid | n | mean skani ANI | skani < 0.95 |
|---|---|---|---|
| [0.95, 0.96) | 85 | 0.9450 | 70 |
| [0.96, 0.97) | 101 | 0.9582 | 4 |
| [0.97, 0.98) | 255 | 0.9703 | 3 |
| [0.98, 0.99) | 369 | 0.9815 | 14 |
| [0.99, 1.01) | 851 | 0.9952 | 0 |

### 2. Chaining inside the largest clusters of species-complex genera

All pairs of up to 60 members (centroid included). `pairs < 95%` counts member pairs that skani puts below 95% although both are within 95% (gmg) of the centroid.

| genus | cluster (centroid) | size in index | members tested | member-centroid < 95% | pairs < 95% | min pair ANI (0 = not reported) | centroid size / median member | species labels |
|---|---|---|---|---|---|---|---|---|
| Aspergillus | GCA_054644345.1_ASM5464434v1 | 518 | 60 | 3/59 | 171/1770 (9.7%) | 0.939 | 1.02 | 5: Aspergillus flavus; Aspergillus oryzae; Aspergillus parasiticus; Aspergillus sojae; Asperg |
| Aspergillus | GCA_058218695.1_ASM5821869v1 | 416 | 60 | 0/59 | 0/1770 (0.0%) | 0.991 | 1.01 | 2: Aspergillus fumigatus; Aspergillus sp. GbtcF1 |
| Aspergillus | GCF_003184595.1_Asplac1 | 101 | 60 | 0/59 | 0/1770 (0.0%) | 0.974 | 1.01 | 3: Aspergillus niger; Aspergillus ochraceus; Aspergillus welwitschiae |
| Aspergillus | GCA_025768885.1_ASM2576888v1 | 83 | 60 | 13/59 | 618/1770 (34.9%) | 0.925 | 0.99 | 12: Alternaria destruens; Aspergillus awamori; Aspergillus costaricensis; Aspergillus luchuens |
| Candida | GCA_042860215.1_ASM4286021v1 | 134 | 60 | 0/59 | 0/1770 (0.0%) | 0.985 | 1.03 | 2: Candida africana; Candida albicans |
| Candida | GCA_032356545.1_MYCT-2423-2021 | 112 | 60 | 0/59 | 0/1770 (0.0%) | 0.977 | 1.01 | 2: Candida sp. LDI48194; Candida tropicalis |
| Candida | GCA_030566835.1_ASM3056683v1 | 15 | 15 | 0/14 | 0/105 (0.0%) | 0.973 | 0.67 | 3: Candida pseudoviswanathii; Candida sp. AC06; Candida viswanathii |
| Candida | GCF_000026945.1_ASM2694v1 | 10 | 10 | 0/9 | 0/45 (0.0%) | 0.991 | 1.02 | 1: Candida dubliniensis |
| Candidozyma | GCA_059629385.1_MRSN455079 | 950 | 60 | 0/59 | 0/1770 (0.0%) | 0.981 | 1.00 | 1: Candidozyma auris |
| Candidozyma | GCF_002926055.2_CanHae_1.0 | 9 | 9 | 0/8 | 0/36 (0.0%) | 0.996 | 1.00 | 1: Candidozyma haemuli |
| Candidozyma | GCA_003708405.2_ASM370840v2 | 6 | 6 | 0/5 | 0/15 (0.0%) | 0.983 | 1.01 | 2: Candidozyma heveicola; Candidozyma molenica |
| Candidozyma | GCF_002926085.2_CanDuoHae_v1.0 | 5 | 5 | 0/4 | 0/10 (0.0%) | 0.999 | 0.99 | 2: Candidozyma duobushaemuli; Candidozyma haemuli |
| Cryptococcus | GCA_056622635.1_ASM5662263v1 | 187 | 60 | 1/59 | 14/1770 (0.8%) | 0.947 | 0.98 | 1: Cryptococcus neoformans |
| Cryptococcus | GCF_000835755.1_ASM83575v2 | 22 | 22 | 13/21 | 135/231 (58.4%) | 0.938 | 1.01 | 6: Cryptococcus bacillisporus; Cryptococcus decagattii; Cryptococcus deuterogattii; Cryptococ |
| Cryptococcus | GCF_000149385.1_ASM14938v1 | 13 | 13 | 0/12 | 0/78 (0.0%) | 0.963 | 0.98 | 2: Cryptococcus deneoformans; Cryptococcus neoformans x Cryptococcus deneoformans |
| Cryptococcus | GCF_001720205.1_Cryp_amyl_CBS6039_ | 5 | 5 | 3/4 | 7/10 (70.0%) | 0.943 | 1.00 | 3: Cryptococcus amylolentus; Cryptococcus floricola; Cryptococcus wingfieldii |
| Fusarium | GCA_040822265.1_ASM4082226v1 | 810 | 60 | 1/59 | 59/1770 (3.3%) | 0.935 | 0.90 | 4: Fusarium foetens; Fusarium odoratissimum; Fusarium oxysporum; Fusarium sp. MPI-SDFR-AT-007 |
| Fusarium | GCA_018345435.1_ASM1834543v1 | 466 | 60 | 2/59 | 116/1770 (6.6%) | 0.941 | 0.99 | 8: Fusarium asiaticum; Fusarium cerealis; Fusarium culmorum; Fusarium graminearum; Fusarium m |
| Fusarium | GCA_022627105.1_ASM2262710v1 | 84 | 60 | 58/59 | 448/1770 (25.3%) | 0.932 | 0.98 | 13: Aspergillus subversicolor; Fusarium annulatum; Fusarium chuoi; Fusarium concentricum; Fusa |
| Fusarium | GCA_030719135.1_ASM3071913v1 | 80 | 60 | 0/59 | 0/1770 (0.0%) | 0.989 | 0.99 | 2: Fusarium poae; Fusarium sp. G73 |
| Penicillium | GCA_023626635.1_ASM2362663v1 | 119 | 60 | 0/59 | 0/1770 (0.0%) | 0.958 | 0.97 | 8: Chaetomiaceae sp. 537B; Penicillium cf. chrysogenum; Penicillium chrysogenum; Penicillium  |
| Penicillium | GCA_052323515.1_ASM5232351v1 | 109 | 60 | 0/59 | 0/1770 (0.0%) | 0.951 | 1.07 | 7: Penicillium biforme; Penicillium camemberti; Penicillium caseifulvum; Penicillium commune; |
| Penicillium | GCA_023141515.1_ASM2314151v1 | 51 | 51 | 0/50 | 0/1275 (0.0%) | 0.956 | 1.05 | 3: Penicillium carneum; Penicillium psychrosexuale; Penicillium roqueforti |
| Penicillium | GCA_001950535.1_JCM_22607_assembly | 36 | 36 | 0/35 | 0/630 (0.0%) | 0.992 | 1.05 | 7: Aspergillus foetidus; Penicillium citrinum; Penicillium hetheringtonii; Penicillium sp. B9 |
| Saccharomyces | GCA_003277715.1_ASM327771v1 | 1302 | 60 | 0/59 | 0/1770 (0.0%) | 0.977 | 1.05 | 1: Saccharomyces cerevisiae |
| Saccharomyces | GCA_002375215.1_ASM237521v1 | 83 | 60 | 4/59 | 277/1770 (15.6%) | 0.000 | 1.05 | 9: Saccharomyces cerevisiae; Saccharomyces cerevisiae x Saccharomyces eubayanus; Saccharomyce |
| Saccharomyces | GCA_030035115.1_ASM3003511v1 | 56 | 56 | 0/55 | 21/1540 (1.4%) | 0.949 | 0.98 | 2: Saccharomyces cariocanus; Saccharomyces paradoxus |
| Saccharomyces | GCA_918298155.1_yHAB447 | 52 | 52 | 5/51 | 186/1326 (14.0%) | 0.938 | 1.00 | 4: Saccharomyces eubayanus x Saccharomyces uvarum; Saccharomyces pastorianus; Saccharomyces u |
| Trichophyton | GCA_056221975.1_ASM5622197v1 | 102 | 60 | 1/59 | 59/1770 (3.3%) | 0.942 | 1.00 | 6: Trichophyton equinum; Trichophyton indotineae; Trichophyton interdigitale; Trichophyton me |
| Trichophyton | GCA_051620515.1_ASM5162051v1 | 47 | 47 | 0/46 | 0/1081 (0.0%) | 0.965 | 1.00 | 3: Trichophyton benhamiae; Trichophyton concentricum; Trichophyton verrucosum |
| Trichophyton | GCA_001651435.1_CMCC_F_T3l | 33 | 33 | 0/32 | 0/528 (0.0%) | 0.990 | 0.99 | 6: Trichophyton kuryangei; Trichophyton rubrum; Trichophyton soudanense; Trichophyton violace |

### 3. Do separate clusters of one genus merge at 95% under skani?

All centroid pairs of the clusters of a genus (up to 150 centroids). gmg's greedy clustering keeps centroids below 95% containment ANI of each other; skani may disagree.

| genus | centroids | pairs | skani ANI >= 95% | >= 97% | size ratio of those pairs (median, larger/smaller n_hashes) | same for all pairs | example (ANI) |
|---|---|---|---|---|---|---|---|
| Aspergillus | 131 | 8515 | 5 | 4 | 1.92 | 1.17 | Aspergillus latus / Aspergillus nidulans (0.999); Aspergillus latus / Aspergillus spinulosporus (0.998) |
| Fusarium | 121 | 7260 | 19 | 16 | 2.28 | 1.19 | Fusarium cf. verticillioides / Fusarium verticillioides (0.988); Fusarium tucumaniae / Fusarium virguliforme ( |
| Saccharomyces | 17 | 136 | 14 | 10 | 3.19 | 1.92 | Saccharomyces kudriavzevii / Saccharomyces kudriavzevii (0.998); Saccharomyces cerevisiae / Saccharomyces sp.  |
| Candida | 17 | 136 | 2 | 2 | 1.70 | 1.13 | Candida tropicalis / Candida tropicalis (0.996); Candida albicans / Candida albicans (0.995) |
| Candidozyma | 11 | 55 | 0 | 0 | - | 1.04 |  |
| Cryptococcus | 14 | 91 | 0 | 0 | - | 1.14 |  |
| Penicillium | 124 | 7626 | 1 | 1 | 4.91 | 1.13 | Penicillium cf. chrysogenum / Penicillium chrysogenum (0.998) |
| Trichophyton | 3 | 3 | 0 | 0 | - | 1.02 |  |

