# F2: GmG marker mode on FusariumID (held-out calibration)

Follow-up 3 of `../F1_fusariumid/README.md`: can gmg identify a fusarioid isolate from marker sequences, the way
FusariumID does, by k-mers instead of alignments? User doc: `docs/marker-mode.md`.

**Index.** One reference per FusariumID strain holding all its loci (tef1, rpb1, rpb2, CaM, tub2, ITS), sketched at
scaled 1 (every k-mer), lineage Ascomycota/Sordariomycetes/Hypocreales/Nectriaceae/genus/species
(`scripts/fusariumid_markers.py prepare`). Full index: 4,851 strains, 676k distinct k-mers, 73 MB, built in about a
minute; held-out index: 3,881 strains.

**Test.** 20% of strains (970, seed 1) removed from the index; each of their marker sequences queried on its own
(3,408 queries; 187 of them from species with no strain left in the index). `run.sbatch` runs gmg with permissive
thresholds and sweeps the species ANI cutoff (`results/calibration.tsv`), then runs the shipped set
`--thresholds markers` (`results/shipped.tsv`, cutoff column irrelevant there: read the 0.98 row).
Run: `REPO=<checkout> W=<build dir> FID=<F1 build>/fid sbatch bench/F2_marker_mode/run.sbatch` (about 1 minute).

## Choosing the species cutoff (pooled over loci)

| species ANI cutoff | species called | precision | species absent from the index named |
|---|---|---|---|
| 0.98 | 0.614 | 0.967 | 0.171 |
| 0.99 | 0.605 | 0.972 | 0.102 |
| **0.995 (shipped)** | **0.588** | **0.974** | **0.048** |
| 0.998 | 0.537 | 0.979 | 0.027 |
| 0.999 | 0.488 | 0.981 | 0.016 |

Precision hardly depends on the cutoff: gmg's vote-consistency gate already withholds the species when the nearest
strains belong to different species (identical or near-identical markers). The cutoff mostly decides what happens to a
species that is not in the index; 0.995 cuts those false names from 17% to 5% for 3 points of called.

## Shipped thresholds (`--thresholds markers`), per locus

| locus | queries | species called | species precision | absent named | FOSC complex calls (precision) | genus called (precision) |
|---|---|---|---|---|---|---|
| tef1 | 930 | 0.725 | 0.977 | 0.109 | 104 (0.981) | 0.982 (1.000) |
| rpb2 | 806 | 0.633 | 0.996 | 0.067 | 87 (1.000) | 0.996 (0.999) |
| rpb1 | 458 | 0.543 | 0.950 | 0 | 36 (1.000) | 0.989 (1.000) |
| CaM | 456 | 0.530 | 0.961 | 0.042 | 63 (0.984) | 0.993 (1.000) |
| tub2 | 331 | 0.594 | 0.973 | 0 | 44 (0.955) | 0.994 (0.994) |
| ITS | 427 | 0.302 | 0.949 | 0 | 23 (1.000) | 0.960 (0.995) |
| all | 3,408 | 0.588 | 0.974 | 0.048 | 357 (0.986) | 0.986 (0.999) |

"species called" is over queries whose species is in the index; FOSC queries that get the complex are counted under the
complex, not as species calls. Compared with F1's alignment nearest neighbour (tef1 79% top-hit correct, ties counted
as not correct), marker mode calls fewer species from one locus but is right 97% of the time when it does, and says
"genus only" or "complex" instead of guessing between sister species.

## Limits

* One locus per query in this test; several loci of one isolate in one FASTA give more species calls (not measured here).
* The held-out split is by strain at random, so close relatives of a query are usually still in the index; the "absent"
  rows are the novelty test (187 queries).
* FusariumID is curated by fusarium.org: the index is not distributed with gmg; users build it from the download.
