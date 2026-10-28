# cluster analysis

23669 genomes, 6104 clusters; singletons: 4155 (68% of clusters, 18% of genomes); largest: 1302

cluster size median / p90 / p99: 1 / 4 / 47

**purity** (genomes whose species label is the most common label in their cluster): 0.848
purity over clusters with >= 2 genomes: 0.815

species with >= 5 genomes: 534; median share of a species in its largest cluster: 1.00; species split over >1 substantial cluster (largest holds < 80%): 93

clusters in which a second species label has >= 5 genomes: 45

| size | centroid | labels |
|---|---|---|
| 810 | GCA_040822265.1_ASM4082226v1 | Fusarium oxysporum (782); Fusarium odoratissimum (8); Fusarium commune (6); Fusarium veterinarium (3) |
| 518 | GCA_054644345.1_ASM5464434v1 | Aspergillus flavus (327); Aspergillus oryzae (128); Aspergillus parasiticus (44); Aspergillus sojae (6) |
| 489 | GCA_020875535.1_ASM2087553v1 | Macrophomina phaseolina (464); Macrophomina pseudophaseolina (12); Macrophomina tecta (9); Macrophomina euphorbiicola (2) |
| 466 | GCA_018345435.1_ASM1834543v1 | Fusarium asiaticum (257); Fusarium graminearum (134); Fusarium culmorum (14); Fusarium boothii (8) |
| 210 | GCA_056882405.1_KAS6096_auto_new | Alternaria alternata (145); Alternaria gaisen (18); Alternaria arborescens (10); Alternaria tenuissima (9) |
| 172 | GCA_023012375.1_ASM2301237v1 | Komagataella phaffii (133); Komagataella pastoris (33); Komagataella kurtzmanii (6) |
| 119 | GCA_023626635.1_ASM2362663v1 | Penicillium chrysogenum (89); Penicillium rubens (14); Penicillium dipodomyicola (3); Penicillium griseoroseum (3) |
| 109 | GCA_052323515.1_ASM5232351v1 | Penicillium commune (65); Penicillium camemberti (15); Penicillium palitans (14); Penicillium fuscoglaucum (6) |
| 104 | GCA_038023885.1_ASM3802388v1 | Colletotrichum siamense (27); Colletotrichum fructicola (26); Colletotrichum gloeosporioides (9); Colletotrichum asianum (4) |
| 102 | GCA_056221975.1_ASM5622197v1 | Trichophyton indotineae (44); Trichophyton mentagrophytes (25); Trichophyton interdigitale (23); Trichophyton tonsurans (5) |
| 101 | GCF_003184595.1_Asplac1 | Aspergillus niger (90); Aspergillus welwitschiae (8); Aspergillus phoenicis (1); Aspergillus awamori (1) |
| 89 | GCA_047371005.1_CABS03 | Colletotrichum lupini (32); Colletotrichum nymphaeae (30); Colletotrichum filicis (8); Colletotrichum acutatum (4) |
| 84 | GCA_022627105.1_ASM2262710v1 | Fusarium proliferatum (32); Fusarium fujikuroi (28); Fusarium annulatum (8); Fusarium oxysporum (3) |
| 83 | GCA_025768885.1_ASM2576888v1 | Aspergillus niger (32); Aspergillus tubingensis (23); Aspergillus luchuensis (14); Aspergillus costaricensis (3) |
| 83 | GCA_002375215.1_ASM237521v1 | Saccharomyces cerevisiae x Saccharomyces eubayanus (37); Saccharomyces pastorianus (16); Saccharomyces eubayanus x Saccharomyces uvarum (9); Saccharomyces cerevisiae x Saccharomyces kudriavzevii (8) |

## leave-one-out, queries whose cluster has >= 2 members (n=19514)

- best hit is in the query's own cluster: 0.982
- species label of the vote equals the query's label: 0.842
- of the queries with a wrong species label, the best hit is still in the query's cluster: 0.963 (n=3081)

## genome-size outliers

69 genomes are > 1.6x the median size (sketch hashes) of their species (species with >= 5 genomes): possible mixed/contaminated assemblies, or true polyploid/hybrid genomes.

| id | species | x median | hashes |
|---|---|---|---|
| GCA_029292925.1_ASM2929292v1 | Ophiocordycipitaceae sp. | 3.31 | 3190 |
| GCA_029292995.1_ASM2929299v1 | Ophiocordycipitaceae sp. | 3.23 | 3113 |
| GCA_029300675.1_ASM2930067v1 | Saccharomycetaceae sp. | 3.21 | 3411 |
| GCA_029303725.1_ASM2930372v1 | Saccharomycetaceae sp. | 3.15 | 3341 |
| GCA_058127425.1_B06 | Buchwaldoboletus xylophilus | 2.86 | 2685 |
| GCA_058127345.1_B25_55_B03 | Buchwaldoboletus xylophilus | 2.8 | 2629 |
| GCA_058129735.1_B07 | Buchwaldoboletus xylophilus | 2.8 | 2636 |
| GCA_058127325.1_B25_55_B04 | Buchwaldoboletus xylophilus | 2.79 | 2626 |
| GCA_058127365.1_B25_55_B05 | Buchwaldoboletus xylophilus | 2.77 | 2603 |
| GCA_058127405.1_B25_55_B01 | Buchwaldoboletus xylophilus | 2.75 | 2583 |
| GCA_000697275.1_RhiMicB9738-1.0 | Rhizopus microsporus | 2.73 | 6274 |
| GCA_054960005.1_ASM5496000v1 | Buchwaldoboletus xylophilus | 2.65 | 2490 |
| GCA_058127385.1_B25_55_B02 | Buchwaldoboletus xylophilus | 2.62 | 2466 |
| GCA_029292965.1_ASM2929296v1 | Ophiocordycipitaceae sp. | 2.36 | 2275 |
| GCA_003314755.1_ASM331475v1 | Pleurotus eryngii | 2.33 | 8911 |
| GCA_029295915.1_ASM2929591v1 | Ascomycota sp. | 2.25 | 4737 |
| GCA_019827475.1_ASM1982747v1 | Penicillium cf. chrysogenum | 2.22 | 6685 |
| GCA_041430195.1_ASM4143019v1 | Agaricus bisporus | 2.18 | 5952 |
| GCA_055773995.1_ASM5577399v1 | Malassezia globosa | 2.12 | 1949 |
| GCA_029294875.1_ASM2929487v1 | Ascomycota sp. | 2.02 | 4248 |
| GCA_003327635.1_ASM332763v1 | Saccharomyces kudriavzevii | 2.01 | 1138 |
| GCF_947243775.1_Skud-IFO1802 | Saccharomyces kudriavzevii | 2.0 | 1129 |
| GCA_918252685.1_IFO1802 | Saccharomyces kudriavzevii | 1.99 | 1125 |
| GCA_918282095.1_yHKS300 | Saccharomyces kudriavzevii | 1.98 | 1116 |
| GCA_980287995.1_PYCC9246 | Zygosaccharomyces siamensis | 1.95 | 1875 |
