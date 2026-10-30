# Real SRA/ENA reads of BFD genomes: depth, paired against single, own genome / species / genus hidden (#42)

100 Illumina WGS runs of fungi that are in the reference, from public archives, classified from the first 5,000 ... 500,000 read pairs with the run's own genome (`self`),
its species or its genus hidden from the reference. Accession list: `data/validation/sra_runs.tsv`; all depth tables: `data/validation/sra_depth_tables.md`;
scripts: `scripts/sra_plan.py`, `sra_harvest.py`, `sra_depth_table.py`, `scripts/hpcc/30_sra_plan.sbatch`, `31_sra_harvest.sbatch`, `32_sra_titration.sbatch`. Per-run results stay
on the HPCC (`GmG_runs/reads_sra`); the raw reads were deleted after scoring.

## Method

* **Which runs.** The ENA Portal API (`read_run`, metadata only, fetched one `first_public` year at a time) lists 236,453 fungal Illumina WGS paired runs. A run is kept when its BioSample is the
  BioSample of a genome of the reference (NCBI `assembly_summary.txt` of the fungal GenBank assemblies, matched by the accession in the genome id) and the species is named the same way by three records:
  the ENA `scientific_name` of the run, the NCBI organism name of the assembly and the reference label; the name is a binomial (not "sp."), the genome is not on the exclusion list
  (`reports/label_audit/reference_exclusions.txt`), the run has at least 600,000 pairs, a mean read length >= 90 and two FASTQ files. That leaves 9,533 eligible runs (827 genera, 16 phyla); one run per genome;
  genera in random order, at most 5 runs per genus and 25 per submitting center, until 100 runs. **Result: 100 runs, 100 genera, 7 phyla (Ascomycota 66, Basidiomycota 24, Mucoromycota 3, Microsporidia 3,
  Mortierellomycota 2, Rozellomycota 1, Neocallimastigomycota 1), 59 submitting centers** (the `center_name` of ENA: 15 runs have none recorded, and many are NCBI submission accessions such as SUB10005639 standing in
  for a lab; JGI 6, University of California Riverside 1, Universidad de Chile 1, RIKEN CLST 2). Read length 90-301 bp (median 150); genome size 2.3-235 Mb (median 36 Mb).
* **Reads.** Only the first 500,000 read pairs of each run were streamed from the ENA file server over HTTPS (`ftp.sra.ebi.ac.uk`; the connection is closed when 2,000,000 lines are decompressed; a dropped
  connection resumes with an HTTP Range request at the next compressed byte). **Route: ENA HTTPS for all 100 runs; the SRA `fastq-dump` fallback was not needed.** **Transferred: 9.57 GB of compressed
  bytes read for 100 of 100 runs (0 failed, 0 short)**, against a limit of 25 GB (about 96 MB per run). A first attempt of the harvester treated an early end of a stream as the end of a file and wrote runs with 1,435 and
  1,994 pairs; it was cancelled after 3 runs (a few hundred MB, deleted, not counted above) and fixed; the figures here come from the second harvest, which accepts a run only when both mates reach 500,000 pairs.
  The per-run byte table was deleted with the raw reads; the total is from the harvest log.
* **Classification.** `gmg run --sample-sheet` with `--exclude-level self|species|genus` and the `ref` column set to the run's own genome id (#85), first N pairs per run (N = 5,000, 15,000, 50,000,
  150,000, 500,000), R1 and R2 merged into one sketch, and R1 only at N = 15,000, 150,000 and 500,000; `known` thresholds, default read settings (minimum abundance 2), index 2026.10.4 (reference
  without the 57 ambiguous-label genomes, with decoys). Binary `gmg 0.0.1` commit `a9acf5cc9` (deployed on the HPCC); **this build predates the `weak_phylum` flag, which is on `main`**. Coverage of an input = bases read /
  the size of the run's own assembly (median 0.04x at 5,000 pairs, 0.12x at 15,000, 0.42x at 50,000, 1.25x at 150,000, 4.2x at 500,000; only genomes up to about 30 Mb exceed 10x, so coverage and clade are confounded in a pooled coverage
  table: the tables below are by N, with the coverage as a median).
* **Truth and what can be found.** Truth = the lineage of the reference genome of the same BioSample (taxdump-normalised). Recall is counted **over attainable inputs only**: with the run's own genome hidden the genus is attainable if
  another genome of that genus exists (70 of the 100 runs) and the species if another genome of the species exists (44); with the species hidden the genus is attainable if another species of the genus exists (61) and the
  species never; with the genus hidden neither. Precision is over the calls made.

## Own genome hidden (another genome of the species or genus may still be in the reference), R1 + R2

| read pairs (median coverage) | genus called | genus precision | genus recall (70 attainable) | genus or tentative genus right | species called | species precision | species recall (44 attainable) | phylum precision |
|---|---|---|---|---|---|---|---|---|
| 5,000 (0.04x) | 0.40 | 0.925 (40) | 0.53 | 0.53 | 0.21 | 0.952 (21) | 0.45 | 1.000 (40) |
| 15,000 (0.12x) | 0.49 | 0.857 (49) | 0.60 | 0.64 | 0.23 | 0.957 (23) | 0.50 | 1.000 (55) |
| 50,000 (0.42x) | 0.50 | 0.840 (50) | 0.60 | 0.69 | 0.24 | 1.000 (24) | 0.55 | 0.956 (68) |
| 150,000 (1.25x) | 0.50 | 0.840 (50) | 0.60 | 0.79 | 0.25 | 0.960 (25) | 0.55 | 0.956 (90) |
| 500,000 (4.2x) | 0.47 | 0.851 (47) | 0.57 | 0.74 | 0.24 | 0.958 (24) | 0.52 | 0.879 (99) |

**Minimum depth.** Recall of 0.8 is not reached at any depth, for genus or for species, in any mode. What depth does decide: **the genus recall is at its plateau (0.60) from 15,000 pairs (0.12x) and the species recall
(0.55) from 50,000 pairs (0.42x)**; more reads change nothing after that, and the genus information (call or tentative genus) keeps growing to 0.79 at 150,000 pairs (1.25x) and then stays. So **about 1x (150,000 pairs) is the depth to
ask for; 0.1-0.4x already gives the plateau of the calls**. The ceiling is not depth.

## Why the recall stops at about 0.6

Of the 30 inputs whose genus is attainable and not called right at 500,000 pairs (own genome hidden), 25 are `low_confidence`, 1 `no_match` and 4 `ok` with another genus. For 19 of the 30 the genus is withheld because the best ANI is below the
genus cutoff (`ani_below_cutoff`): the nearest remaining genome is another species of the genus at ANI 0.78-0.89 (for example *Epibryon bryophilum* 0.82, *Paxillus involutus* 0.89, *Hermanozyma aechmeae* 0.89), and
the top hit's genus is the right one for 15 of the 30 and the tentative genus for 12. 4 are withheld for `low_agreement` of the near-best hits (one is *Glaciozyma*: ANI 0.9989 to *G. watsonii* but the near-best hits disagree on the family,
0.69 < 0.90: conflicting reference labels). That is the design of the thresholds (calls only above an ANI that has been calibrated for precision), not a lack of reads. A "genus" is attainable only in principle when the congeners
left in the reference are other species; with the species also hidden the genus recall is 0.33-0.36 (table below).

## Species hidden and genus hidden (R1 + R2)

| hidden | read pairs | genus called | genus precision | genus recall (attainable) | genus or tentative genus right | species called (must be wrong) | phylum precision |
|---|---|---|---|---|---|---|---|
| species | 15,000 | 0.34 | 0.647 (34) | 0.36 (of 61) | 0.44 | 0.04 (4) | 1.000 (43) |
| species | 150,000 | 0.33 | 0.636 (33) | 0.34 | 0.62 | 0.04 (4) | 0.943 (88) |
| species | 500,000 | 0.30 | 0.667 (30) | 0.33 | 0.59 | 0.04 (4) | 0.827 (98) |
| genus | 15,000 | 0.12 | 0.000 (12) | (none attainable) | - | 0.02 (2) | 1.000 (19) |
| genus | 150,000 | 0.13 | 0.000 (13) | - | - | 0.02 (2) | 0.819 (83) |
| genus | 500,000 | 0.11 | 0.000 (11) | - | - | 0.02 (2) | 0.745 (94) |

With the species hidden, 4% of the inputs get a species call (necessarily wrong) and about a third get a genus call at precision 0.64-0.67 (a congener or a conflicting label); with the genus hidden 11-13% get a genus call that must be wrong
(the same 12% as for assemblies in the truth-set report and 11.6% for RNA-seq).

## Paired against single (own genome hidden)

| read pairs | R1 + R2: genus recall / genus or tentative / species recall / phylum called | R1 only: same |
|---|---|---|
| 15,000 | 0.60 / 0.64 / 0.50 / 0.55 | 0.51 / 0.51 / 0.55 / 0.41 |
| 150,000 | 0.60 / 0.79 / 0.55 / 0.90 | 0.60 / 0.73 / 0.57 / 0.72 |
| 500,000 | 0.57 / 0.74 / 0.52 / 0.99 | 0.57 / 0.74 / 0.55 / 0.94 |

(phylum called = the share of the 100 runs with a phylum call; the species-hidden mode gives the same picture: genus recall 0.36 / 0.34 / 0.33 for pairs and 0.23 / 0.33 / 0.33 for R1 only.) Genus and species
recall are the same once the plateau is reached; what the second read buys is **speed to the plateau**: at 15,000 reads the pair has 0.60 genus recall where R1 alone has 0.51, and at 150,000 the phylum call exists for 90 runs against 72.
A pair is worth about twice the bases of one file, so the single file needs about twice the reads for the same result, as in the ZyGoLife titration.

## By clade and lab (own genome hidden, 500,000 pairs)

| phylum | runs | genus attainable | genus called | genus precision | genus recall (attainable) | phylum precision |
|---|---|---|---|---|---|---|
| Ascomycota | 66 | 49 | 0.50 | 0.909 (33) | 0.61 | 0.908 (65) |
| Basidiomycota | 24 | 14 | 0.46 | 0.727 (11) | 0.57 | 0.958 (24) |
| Mucoromycota | 3 | 2 | 0.67 | 1.000 (2) | 1.00 | 0.667 (3) |
| Microsporidia | 3 | 3 | 0.00 | - | 0.00 | 0.667 (3) |
| Mortierellomycota | 2 | 2 | 0.00 | - | 0.00 | 0.000 (2) |
| Neocallimastigomycota, Rozellomycota | 1 + 1 | 0 | - | - | - | - |

Per genus there is one run, so a per-genus table is the accession list. Per lab (centers with at least 3 runs; the lab is a submitter accession, not a person): `?` (15 runs) genus recall 0.60 (3 of 5), SUB10005639 (8) 0.25, JGI (6) 0.00 of 3
attainable, SUB3485651 (5) 0.33, SUB10689706 (5) 0.50, SUB14100067 (5) 0.00 of 1 attainable: counts of 1-8, so the spread (0.00-0.60) is noise plus the thin clades each submitter sequences, and no lab effect can be claimed. Microsporidia and Mortierellomycota are not called to genus at 500,000 pairs
(3 + 2 runs; the other two thin phyla have no attainable genus); Ascomycota and Basidiomycota behave alike.

## Confident calls

* **Wrong phylum calls are all `low_confidence`:** 12 of 99 phylum calls at 500,000 pairs (own genome hidden) and 17 of 98 (species hidden), every one `low_confidence`, none among the 47 (30) `ok` calls. Most are inputs from clades the reference holds
  few genomes of whose best hit has ANI 0.75-0.83 to a few large repetitive genomes (*Gigaspora*, *Rhizophagus*, *Puccinia*, *Massospora*); the decoy rerun (#89) keeps them from becoming a non-fungal phylum. The fall of the
  phylum precision with depth (1.000 at 5,000 pairs, 0.879 at 500,000) is the extra `low_confidence` calls that deeper reads allow, not worse calls.
* **Wrong genus calls among the 47 `ok` calls: 7**, all taxonomy or label questions on inspection (my reading, not an authority): 3 are genera absent from the rest of the reference, named after the nearest genus (*Gastroboletus* -> *Boletus* ANI 0.91, *Anaeromyces* ->
  *Neocallimastix* 0.92, *Sistotremastrum* -> *Sertulicium* 0.97); *Talaromyces purpureogenus* -> *Penicillium* 0.987 and *Lecanicillium aphanocladii* -> *Aphanocladium album* 0.990 are generic synonymy; *Xerocomellus* -> *Boletus* 0.969 (Boletaceae).
  **One is a probable mislabel: ERR14883901, registered as *Coccidioides immitis* (5 other *C. immitis* genomes remain in the reference), has *Ophidiomyces ophidiicola* as best hit at ANI 0.947**: either the run or the reference genome of that BioSample is
  not *C. immitis* (for #14).
* The species calls that are made are precise: 0.95-1.00 at every depth (24-25 calls at 50,000 pairs and more).

## Comparison with the other real-reads sets

| set | what | genus | species |
|---|---|---|---|
| this report, own genome hidden (100 runs, 100 genera) | genus attainable for 70 | called 0.47-0.50, precision 0.84-0.86, recall 0.57-0.60 | called 0.24, precision 0.96, recall 0.52-0.55 (44 attainable) |
| ZyGoLife (119 early-diverging fungi, `reads_titration.md`) | genus in the reference for 98 | called 0.78 from 0.5x; precision of `ok` calls 0.72-0.76 against lab labels (7 of 8 checked misses are the label) | not scored |
| *Coccidioides* (558 sets, `cocci_reads_*.md`) | species and its relatives are in the reference | called >= 99.6%, precision 0.991-0.993 | called 95-96%, precision 0.994-0.996 |
| RNA-seq (345 sets, `reads_rna_exclusion.md`) | species hidden | called 0.47, precision 0.875 | called 0.10, precision 0 (must be wrong) |

The three real-DNA sets are consistent once the reference content is taken into account: reads call a genus and species at assembly-like precision (0.84-1.00) when a close relative is in the reference (*Coccidioides*: hundreds of conspecific genomes; ZyGoLife with
the genus present: called 0.87), and the calls are withheld, not wrong, when the nearest genome is another species at ANI < 0.9 (here). The precision in this set is higher than the ZyGoLife one probably because the truth is the reference label of the same BioSample, not a lab
name (not tested). The genus calls stop gaining from more reads at about 0.5-1x in the ZyGoLife and in this set; the *Coccidioides* depth series shows species calls still rising to 2M reads.

## Acceptance of #42

At least 5 genera from at least 3 labs: **100 genera from 59 submitting centers, 7 phyla**. Depth curves with the minimum depth for genus and species: given above; recall 0.8 is not reached in any mode, the plateau is at 0.12x (genus) and 0.42x (species), and the
genus information reaches 0.79 at about 1x. Paired against single: given above. The reads are real public data (ENA), 100 of 100 runs, 9.57 GB.

## Not covered

Read length < 100 bp or single-end libraries (the selection is paired Illumina >= 90 bp), long reads (#44), runs of fewer than 600,000 pairs, coverage above about 10x for genomes over 30 Mb (500,000 pairs is the cap),
and a curated truth beyond the reference label of the BioSample's own assembly.
