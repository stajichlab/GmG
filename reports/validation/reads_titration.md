# DNA reads: coverage titration on real paired Illumina reads (#42)

Reads of ZyGoLife strains (`/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Assembly/input`, read-only: 975 paired NovaSeq 2 x 150 bp samples, 76 genera, mostly
Mucoromycota and other early-diverging fungi, which are thin in the reference). Two strains per genus (`scripts/make_reads_sheet.py`): 120 read sets, 119 scored; truth = the NCBI
lineage of the genus in the file name (`scripts/name_truth.py`, old genera resolved through the taxdump synonyms). **The truth is the lab's label, not a curated identification**:
see "Labels" below. Truth phyla: Mucoromycota 67, Chytridiomycota 21, Kickxellomycota 15, Mortierellomycota 8, Blastocladiomycota 3, Basidiobolomycota 2, Entomophthoromycota 2, Ascomycota 1.
The truth genus is in the reference for 98 of the 119 (`genus_in_ref`), not in it for 21. Species is not scored (lab species names are not curated).

Each strain is classified from the first N read pairs (R1 and R2 merged into one sketch by the sample sheet), N = 3,000 ... 1,500,000 (about 0.02x ... 10x of a 45 Mb genome),
index 2026.10.3, `known` thresholds (`scripts/hpcc/25_reads_titration.sbatch`, scoring with `gmg validate`). Per-input results on the HPCC (`GmG_runs/reads/titration`).

Cells: *called share / precision (calls made)*.

## All inputs

| read pairs | n | phylum | class | order | family | genus |
|---|---|---|---|---|---|---|
| 3,000 | 119 | 0.39 / 0.913 (46) | 0.39 / 0.913 (46) | 0.39 / 0.913 (46) | 0.39 / 0.826 (46) | 0.39 / 0.761 (46) |
| 15,000 | 119 | 0.75 / 0.910 (89) | 0.75 / 0.910 (89) | 0.75 / 0.898 (88) | 0.75 / 0.818 (88) | 0.73 / 0.736 (87) |
| 75,000 | 119 | 0.87 / 0.913 (103) | 0.87 / 0.913 (103) | 0.86 / 0.902 (102) | 0.83 / 0.775 (98) | 0.78 / 0.720 (93) |
| 300,000 | 119 | 0.93 / 0.919 (111) | 0.93 / 0.919 (111) | 0.92 / 0.917 (109) | 0.91 / 0.785 (107) | 0.78 / 0.742 (93) |
| 1,500,000 | 119 | 0.99 / 0.898 (118) | 0.99 / 0.898 (118) | 0.95 / 0.902 (112) | 0.91 / 0.794 (107) | 0.78 / 0.753 (93) |
| 75,000, R1 only | 119 | 0.78 / 0.925 (93) | 0.78 / 0.925 (93) | 0.78 / 0.913 (92) | 0.77 / 0.802 (91) | 0.74 / 0.739 (88) |

## Confident calls only (status `ok`) and the rest (`low_confidence`)

| read pairs | `ok` n | phylum | genus | `low_confidence` n | phylum | family |
|---|---|---|---|---|---|---|
| 3,000 | 46 | 0.913 | 0.761 | 0 | - | - |
| 15,000 | 87 | 0.908 | 0.736 | 2 | 1.000 | 0.500 |
| 75,000 | 93 | 0.925 | 0.720 | 10 | 0.800 | 0.333 |
| 300,000 | 93 | 0.946 | 0.742 | 18 | 0.778 | 0.533 |
| 1,500,000 | 93 | 0.946 | 0.753 | 25 | 0.720 | 0.600 |

(Precision; every `ok` input has a genus call, no `low_confidence` input does. At 3,000 pairs 74 of the 120 classified inputs are `no_match`: a median of 7 hashes; at 75,000 pairs 16; at 1.5M 1.)

## What it shows

* **Depth gets an input to a call, and it stops mattering early.** The share of inputs with a phylum call rises 0.39, 0.75, 0.87, 0.93, 0.99 from 3k to 1.5M pairs; the share with a genus call
  stops at 0.78 (93 inputs) from 75,000 pairs (0.5x) on, and so does the number of `ok` inputs (93). The remaining 26 inputs never reach a genus call at any depth: their genus is absent or
  thin in the reference (the `genus_in_ref = False` half of the table below), so more reads cannot help; a richer reference can.
* **A merged pair is worth about twice the reads of one file**: at 75,000 pairs R1 alone gives phylum called 0.78 (0.87 with the pair) and genus called 0.74 (0.78), with the same precision.
* **Precision is flat in depth, so it is not a depth effect.** The genus precision of the confident calls is 0.72-0.76 from 3k to 1.5M pairs, the phylum precision 0.91-0.95. Reads and
  assemblies of the same strains give the same answer (next section), so the difference from the label is not noise of the reads.
* **Genus not in the reference** (21 inputs): at 1.5M pairs the phylum is right for 0.762 of the calls (all 21 called) and a genus is called for 8 (0 correct, as it must be: the genus is not there);
  with the genus in the reference (98 inputs) phylum 0.928, genus called 0.87 and right for 0.824.

## Labels: reads against the assembly of the same strain

The assemblies of these strains were classified in the external validation (#40, `GmG_runs/external/zygolife`). For the 46 read sets whose strain has an assembly (the strain id is the end of the file name;
`scripts/reads_vs_assembly.py`), at 1.5M pairs: **the genus call of the reads equals the genus call of the assembly for 27 of 28 inputs with both called (0.964) and the phylum call for 42 of 46 (0.913)**.
Where the read call differs from the lab label (8 genus calls), **the assembly makes the same call as the reads for 7**:

| lab label | reads | assembly |
|---|---|---|
| *Absidia* | *Lichtheimia* | *Lichtheimia* |
| *Chaetocladium* | *Cunninghamella* | *Cunninghamella* |
| *Circinomucor* | *Mucor* | *Mucor* |
| *Conidiobolus* | *Gilbertella* | *Gilbertella* |
| *Conidiobolus* | *Thalia* (a decoy, see below) | *Thalia* |
| *Gongronella* | *Halteromyces* | *Halteromyces* |
| *Rhizomucor* | *Phycomyces* | *Phycomyces* |
| *Syncephalastrum* | *Cunninghamella* | *Syncephalastrum* |

So most of the "wrong" genus calls are the label or the reference being the odd one out (modern genus names such as *Absidia* -> *Lichtheimia*, mislabelled strains or reference genomes, poorly sampled
clades), not a failure of the reads method. This was checked for the 46 matched strains only; the other 73 were not reviewed one by one, so the 0.75 is a lower bound on the precision of
the confident genus calls and must not be quoted as the accuracy of GmG on reads. Cases worth a label check (#14): the *Mucor* read set with a *Irpex lacteus* (Basidiomycota) best hit at ANI 0.94,
*Nowakowskiella* (Chytridiomycota) with a *Gonapodya* (Monoblepharomycota) best hit at ANI 0.99, and *Paraphysoderma* with *Boothiomyces* at 0.98: ANI that high between a "wrong phylum" pair says the lab or reference label is wrong.

## A problem for the decoy gate: *Thalia democratica* and the Entomophthoromycota

One `ok` read set (*Conidiobolus antarcticus* ARSEF 691, Entomophthoromycota) gets the phylum Chordata and the flag `not_fungal`: its best hit is the decoy *Thalia democratica* (a salp,
GCF_965202585.1) at ANI 0.913 (182 shared hashes, 15% of the query). It is not an accident of the reads. In the assembly validation (#40) 9 of the 552 ZyGoLife assemblies have that decoy as
the best hit: 7 *Conidiobolus*/*Neoconidiobolus* (ANI 0.85-0.91) and 2 *Mucor hiemalis* (0.83); two of them cross the gate and are flagged `not_fungal` (ARSEF 6913 at ANI 0.9001, ARSEF 79 at 0.9055; shared
hashes 825 and 982, 11-12% of the assembly). Chance cannot give a salp genome 10% of the k-mers of a fungus (a 4^21 k-mer space against a few hundred Mb), so the decoy shares
real sequence with these fungi: probably fungal contaminant contigs in the *Thalia* assembly, or a conserved repeat. The consequence is that the `not_fungal` gate (ANI >= 0.90), set from test sets that had
no Entomophthoromycota, gives a **false `not_fungal` on a real fungus for the thin clade Entomophthoromycota** (3 inputs here: 2 assemblies and 1 read set, all *Conidiobolus*-type). **Test** (`scripts/hpcc/28_thalia_check.sbatch`: the 9 assemblies and the reads set classified with and without that decoy, `--exclude-ids`, index 2026.10.3): hiding the decoy removes
the false `not_fungal` and the Chordata phylum for all three inputs (ARSEF 6913, ARSEF 79 and the reads set); they become `low_confidence` calls with a fungal best hit at ANI 0.88-0.90
(*Neocallimastix* in the reference, the nearest genus it has). The other six inputs do not change. Note that the *Thalia* genome was closer to these fungi (ANI 0.90-0.91) than any fungal genome in the reference
(0.88-0.90): the reference has no close Entomophthoromycota. For the parent project: remove GCF_965202585.1 from the decoy panel (the panel has other animals; one salp is not needed), or check that assembly
for fungal contigs, before the next release; this looks like a problem of that decoy, not of the 0.90 threshold: raising the gate to 0.92 would hide these three, but only by luck, because the decoy is closer than any fungal genome in the reference.
Whether the *Thalia* assembly contains fungal contigs was not examined.

## Not covered

Real SRA reads of other clades (the ZyGoLife and *Coccidioides* sets are all from this lab; the independent assemblies are in #40), single-end short reads below 100 bp, long reads (#44), and a species truth
(lab species names are not curated). The *Coccidioides* depth series (558 sets, species 0.996) is in `cocci_reads_*.md`.
