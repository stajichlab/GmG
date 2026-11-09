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

RESULTS_PLACEHOLDER
