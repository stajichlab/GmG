# F1: how much information separates fusarioid species, and how far can GmG go? (FusariumID vs genome sketches)

Question (2026-11-08): FusariumID (https://www.fusarium.org/page/ID) identifies fusarioid isolates from a few marker loci.
How much sequence does it take to separate species, and to what extent could GmG (whole-genome k-mer sketches) do the
same job?

Data: the FusariumID download `sequencesindatabase23062026.zip` (one FASTA, `fusarioidID23062026.fas`: 18,844
sequences, 4,855 strains, 709 species after name normalisation, 20 genera; loci tef1 4,658, rpb2 4,063, CaM 2,443, rpb1
2,280, ITS 2,118, tub2 1,688, LSU 1,031, minor acl1/his3/IGS/actin). Genomes: the 2,111 BFD genomes of fusarioid
genera (2,099 *Fusarium* as NCBI files them, which includes *Neocosmospora*; one file is empty), 301 NCBI species names.

## Steps (all on the HPCC; build dir `~/projects/GmG/bench/F1_fusariumid_build/`, purgeable)

| step | script | what |
|---|---|---|
| 0 | `parse_fusariumid.py` | headers -> `records.tsv` + `loci/<locus>.fa`; species names normalised (synonyms -> accepted name, `sp.` taxa kept whole, annotations dropped) |
| 1 | `run_markers.sbatch` -> `marker_resolution.py` | all-vs-all blastn per locus; leave-one-out nearest-neighbour species call (other strains only), single loci and multilocus sets; differences to the nearest other species; the same with k-mers (k=21, FracMinHash at scaled 1-100) |
| 2 | `run_genomes.sbatch` (`make_lists.sh`, `genome_ani.py`) | `gmg sketch` k21 s1000; all-pairs sketch ANI; `gmg bench` loo/los against the production indexes at scaled 1000, 2000, 10000, 20000, 50000, with the full BFD reference set (`all`) and fusarioid-only (`fus`) |
| 3 | `run_bridge.sbatch` -> `genome_markers.py` | megablast of the tef1/rpb1/rpb2/CaM/tub2 FusariumID sequences against every genome (~5 s/genome); per-locus and multilocus FusariumID species call per genome |
| 4 | `genome_resolution.py` | NCBI vs FusariumID labels; gmg accuracy against both; genome ANI gaps per species complex |

Run: `REPO=<checkout> W=<build> FAS=<fas> sbatch run_markers.sbatch`, then `run_genomes.sbatch` and `run_bridge.sbatch`
(same `REPO`/`W`), then `python3 genome_resolution.py $W $W/genome_results`. Tests: `pixi run pytest tests/python/test_f1_fusariumid.py`.
Result tables copied here: `results/`.

## Answer in short

* **Markers: about 3 loci / ~4 kb are needed, and species are separated by very little sequence.** Leave-one-out
  nearest-strain accuracy is 79% with tef1 alone, 92.5% tef1+rpb2, 96% tef1+rpb1+rpb2, 97% with four loci. Half of all
  species differ from their closest relative by <= 5 differences in tef1 (~640 bp); 21% are identical to another species
  in tef1. With tef1+rpb1+rpb2 the median is 25 differences (~0.6%), and 6% of species are still within 2.
* **Genomes: GmG has far more information than it needs; the limit is the labels, not the sketch.** gmg leave-one-out
  species accuracy is the same from scaled 1000 (44k hashes per genome) to scaled 50000 (884 hashes): 91.5% against NCBI
  names, 94.5% against FusariumID names. At gmg's species gate (top ANI >= 0.9938) 83% of queries are called and
  96.7% of those agree with FusariumID.
* **Labels are the bottleneck.** In-silico FusariumID barcoding of the 2,110 BFD fusarioid genomes gives a confident
  call (<= 1% from a FusariumID strain over tef1/rpb1/rpb2/CaM/tub2) for 2,059; only 897 agree with the NCBI name.
  Most of the rest are NCBI "*F. oxysporum*" split into FOSC species, *F. solani* -> *Neocosmospora*, *F. proliferatum*
  -> *F. annulatum*.
* **FOSC is where both approaches run out.** Genomes of different FusariumID FOSC species reach 99.6-99.99% ANI; 38% of
  FOSC genomes have a nearest other-species genome above gmg's species gate, and only 16% of FOSC species show a genome
  barcode gap. Markers are not much better there (tef1 44% correct, 3 loci 90%). Outside FOSC (and a few
  *Neocosmospora*) no genome's nearest other-species genome passes the species gate.
* **Marker sequences themselves are invisible to GmG's genome indexes**: at scaled >= 1000 a 640 bp tef1 amplicon keeps
  ~0.6 hashes. With all k-mers (scaled 1) k-mer nearest neighbour on tef1 is 92% correct, at scaled 10 (63 hashes) 91%,
  at scaled 100 (8 hashes) 71%. A FusariumID-style "marker mode" would need its own all-k-mer index (17 Mb of sequence).

## Markers (FusariumID only)

Leave-one-out, query vs other strains; answerable = species with >= 2 strains having the loci (`results/markers_*.tsv`).
"ambiguous" = the nearest strains include another species at exactly the same distance (mostly identical sequences).

| loci | strains answerable | correct | ambiguous | wrong | median diffs to nearest other species | species <= 2 diffs |
|---|---|---|---|---|---|---|
| tef1 | 4,466 | 0.795 | 783 | 134 | 5 | 33% |
| rpb2 | 3,864 | 0.671 | 1,154 | 118 | 3 | 49% |
| rpb1 | 2,139 | 0.723 | 394 | 199 | 3 | 48% |
| CaM | 2,289 | 0.590 | 873 | 66 | 2 | 65% |
| tub2 | 1,580 | 0.665 | 457 | 73 | 2 | 60% |
| ITS | 1,924 | 0.422 | 841 | 272 | 0 | 75% |
| LSU | 859 | 0.348 | 471 | 89 | 0 | 84% |
| tef1+rpb2 | 3,774 | 0.925 | 192 | 91 | 11 | 14% |
| tef1+rpb1 | 2,108 | 0.945 | 28 | 88 | 16 | 7% |
| tef1+rpb1+rpb2 | 2,014 | 0.959 | 10 | 72 | 23 | 3% |
| tef1+rpb2+CaM+tub2 | 973 | 0.971 | 5 | 22 | 18 | 3% |
| tef1+rpb1+rpb2+CaM+tub2 | 666 | 0.963 | 0 | 23 | 29 | 0.3% |

Multilocus rows use the strains that have every locus of the set, so n shrinks and the sets are not the same strains.
By complex (tef1+rpb1+rpb2): FSAMSC 98%, FFSC 97%, FIESC 97%, FTSC 95%, *Neocosmospora* 94%, FOSC 90%; tef1 alone in
FOSC 44% (86% of FOSC strains are within 2 tef1 differences of another species). The table's "median diffs" is over
query strains; over species (`markers_nearest_species.tsv`) it is 5 for tef1 and 25 for tef1+rpb1+rpb2.

k-mer nearest neighbour (k=21, Jaccard; `markers_kmer.tsv`), accuracy by FracMinHash scaled (mean hashes per sequence):

| locus | s=1 | s=3 | s=10 | s=30 | s=100 |
|---|---|---|---|---|---|
| tef1 | 0.918 (621) | 0.917 (199) | 0.910 (63) | 0.882 (22) | 0.715 (8) |
| rpb2 | 0.872 (1363) | 0.872 | 0.864 (137) | 0.825 | 0.720 (13) |
| rpb1 | 0.862 (1389) | 0.861 | 0.856 (143) | 0.828 | 0.675 (12) |
| CaM | 0.785 | 0.770 | 0.758 | 0.680 | 0.354 (4) |
| tub2 | 0.805 | 0.803 | 0.781 | 0.681 | 0.481 (8) |
| ITS | 0.656 | 0.633 | 0.605 | 0.533 | 0.278 (7) |

k-mers score higher than the alignment nearest neighbour on single loci because Jaccard breaks the exact-distance ties
that alignment reports as ambiguous (it also rewards length/region overlap, so it is not a cleaner measure of
divergence); the order of loci is the same.

## Genomes (GmG)

FusariumID labels per genome (`genomes_labels.tsv`, `genomes_label_summary.tsv`): 2,059 of 2,110 called within 1%
(1,921 within 0.2%); 143 are ties between FusariumID species (mostly FOSC). Agreement with NCBI: same 897, differs 1,123,
NCBI name `cf.`/`sp.` 38, no FusariumID call 52. Top relabelings: *F. oxysporum* -> *F. cugenangense* 138,
*F. fabacearum* 121, *F. fabacearum/F. wimaladesilvae* (tie) 93, *F. nirenbergiae* 81, *F. odoratissimum* 42;
*F. proliferatum* -> *F. annulatum* 37; *F. solani*/*cf. solani*/*falciforme* -> *Neocosmospora* spp.

gmg bench, fusarioid-only references (identical with the full BFD reference set; `genomes_bench_summary.tsv`):

| scaled | hashes/genome | loo top = NCBI species | loo top = FusariumID species | species gate called | precision vs FusariumID |
|---|---|---|---|---|---|
| 1000 | 44,491 | 0.915 | 0.945 | 0.830 | 0.967 |
| 10000 | 4,514 | 0.913 | 0.943 | 0.832 | 0.966 |
| 50000 | 884 | 0.916 | 0.946 | 0.836 | 0.967 |

Per complex at scaled 10000 (`genomes_complex_bench.tsv`; FusariumID labels; "los above gate" uses bench's NCBI-label
leave-species-out and is inflated by NCBI mislabels, see the ANI table for the FusariumID version):

| complex | genomes | top correct | gate called | gate precision |
|---|---|---|---|---|
| FOSC | 801 | 0.889 | 0.863 | 0.923 |
| FSAMSC | 578 | 0.997 | 0.971 | 1.000 |
| FFSC | 240 | 0.967 | 0.779 | 1.000 |
| *Neocosmospora* | 162 | 0.938 | 0.617 | 0.980 |
| FIESC | 49 | 1.000 | 0.674 | 1.000 |

Genome ANI gap (all-pairs sketch ANI at scaled 1000, max-containment; `genomes_ani_gap.tsv`; genomes of labels with >= 2 genomes):

| labels | genomes | species | nearest genome same species | species with a gap | other-species ANI median / q90 / max | other species above 0.9938 |
|---|---|---|---|---|---|---|
| NCBI | 1,926 | 116 | 0.910 | 0.526 | 0.988 / 0.996 / 1.000 | 14.9% |
| FusariumID | 1,926 | 148 | 0.942 | 0.669 | 0.989 / 0.995 / 1.000 | 16.4% |
| FOSC | 801 | 38 | 0.886 | 0.158 | 0.993 / 0.997 / 1.000 | 38.5% |
| FSAMSC | 578 | 27 | 0.997 | 0.926 | 0.988 / 0.992 / 0.993 | 0 |
| FFSC | 240 | 27 | 0.975 | 0.852 | 0.984 / 0.988 / 0.994 | 0 |
| *Neocosmospora* | 163 | 22 | 0.939 | 0.682 | 0.984 / 0.992 / 1.000 | 4.3% |

## What GmG could do (proposals, not built)

1. **Name Fusarium genomes with FusariumID-consistent labels.** The sketch side is ready (even the 50000-scaled web tier
   resolves these species); the reference labels are not. The genome -> FusariumID table here is a candidate
   relabelling; per the working rules it belongs in PopSeq (metadata source of truth), not in GmG.
2. **Call FOSC at complex level.** Report "F. oxysporum species complex" plus the nearest FusariumID species and its
   ANI, instead of a species name, when the top hit is in FOSC. A species-complex rank (or a per-genus species gate)
   would also fit other genera with the same problem.
3. **Marker mode.** A separate all-k-mer index of the FusariumID loci (scaled 1, ~17 Mb of sequence) would let GmG
   take a Sanger tef1/rpb2 read; tef1 nearest neighbour is ~92% at scaled 1-3 and still 91% at scaled 10.
   Alignment-free calls would be weaker than FusariumID's curated multilocus pipeline for close species (2-5 differences).
4. **In-silico barcoding inside GmG.** Report a genome's FusariumID nearest strain per locus (what `genome_markers.py`
   does with blast in ~5 s/genome) next to the k-mer call; disagreements flag mislabels or FOSC ties.

## Caveats

* FusariumID taxonomy is the yardstick; the genome calls inherit its species concepts (and FOSC's are debated).
* Marker distance = (mismatches + gap openings) / overlap of the best HSP; no coverage floor, because some records join
  two rpb2 amplicon regions (a 50% floor first dropped their conspecific pairs; fixed in `marker_resolution.py`).
  The genome bridge keeps a 60% query-coverage floor, so two-region rpb2 records do not contribute there.
* bench writes the raw species vote; the species-gate numbers apply only gmg's ANI gate (0.9938), not its other gates.
* One BFD fusarioid genome file is empty (`gmg sketch` FAIL 1).
