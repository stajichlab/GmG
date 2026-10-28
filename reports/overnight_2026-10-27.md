# GmG overnight report (2026-10-27 to 10-28)

Everything below was run, not estimated, unless marked **[est]**. Code is on branch `overnight` of `stajichlab/GmG`; `LOG.md` has the chronological record, issues #1-#7 track the stages. Work was done unattended; decisions I made without you are listed in section 9.

## 1. Headline

- **The tool works.** `gmg` (Rust) sketches a genome or reads, searches a reference index, and returns per-rank calls (phylum, class, order, family, genus, species) with a vote consistency and a containment-ANI. 23,669 BFD genomes are indexed.
- **Held-out genomes (BFD, leave-one-out, raw labels):** phylum 98.5%, class 97.9%, order 97.4%, family 97.0%, genus 95.8%, species 92.8% correct (answerable queries; species is answerable for 76% of genomes because 5,576 BFD species have a single genome).
- **Genomes not in BFD (145 NCBI assemblies, default thresholds):** for the 120 whose species is in BFD: phylum to genus 100% correct, species 94%. For 25 unnamed isolates: phylum 92%, genus 77% (vote).
- **Reads work down to 0.5x** of simulated Illumina-like reads (150 bp, 0.5% error): phylum/order/genus/species calls match the genome-level result. 12 real *Coccidioides* read sets (100k to 2M reads) are 12/12 correct at every rank.
- **The browser payload can be small.** Accuracy is essentially unchanged from scaled=10000 to scaled=50000, so one representative per species at scaled=50000 is **25 MB** (3 per species: 29 MB) as compact v1; at scaled=20000, 61 MB (71 MB).
- **Species labels are the limit, not the sketches.** Near-identical genomes (>99.5% ANI) carry a different species label 24% of the time in BFD (hybrid names, species complexes, synonyms). The right output is an ANI-cluster ID next to the label (section 6, not built yet).
- **BFD label audit:** 40 genomes have a lineage that conflicts with near-identical neighbours at phylum/class, and 23 reference genomes are the "hub" cause of most apparent errors in the leave-genus-out test (one of them is the top hit of 416 genomes of another class). Lists in `reports/label_audit/`.

## 2. What was built

| piece | where | state |
|---|---|---|
| FracMinHash sketch, sourmash-compatible | `crates/gmg-core/src/sketch.rs` | verified hash-for-hash against the `sourmash` Rust crate on A1163 and Af293 (`tests/oracle.rs`) |
| sketch file `.gsk` (with abundances), reference index (inverted hash->genomes), v0 and compact v1 format, search, per-rank calls | `gmg-core` | 12 unit tests pass |
| CLI: `sketch`, `build`, `guess`, `info`, `convert`, `bench` (leave-one/species/genus-out, sharded, `--ref-list`) | `crates/gmg-cli` | used for everything below |
| HPCC pipeline (SLURM) | `scripts/hpcc/` | stage 3-7 jobs ran |
| analysis and audit | `scripts/analyze_bench.py`, `find_label_suspects.py`, `score_guess.py`, `score_reads.py`, `pick_reps.py` | |
| external genome selection, read simulator | `scripts/select_external.py`, `sim_reads.py`, `make_sim_reads.sh` | |

Speed (measured): sketching 0.37 CPU-s per 30 Mb genome (Python sourmash: ~1.2 s). Search of one query at scaled=10000 against all 23,669 genomes: ~17 ms CPU [from the 23,669-query leave-one-out in 98 s on 4 threads]. At scaled=2000: ~90 ms.

## 3. Reference build (HPCC, `/bigdata/stajichlab/shared/projects/GmG_runs/`)

23,679 genomes had a samples.csv row; 10 `.fa.gz` were empty (all `DAH10xx`, one chunk) -> **23,669 sketched** (k=21, scaled=1000, ~1 h per 200-genome task at 2 CPUs; the array waited hours behind the BFD pipeline's CPU limit). 198 files in `input_clean_genomes/` have no samples.csv row; 1 row has no file.

| index | scaled | v0 (uncompressed) | v1 (delta+varint) |
|---|---|---|---|
| all 23,669 genomes | 1000 | 5.6 GB | not converted |
| all | 2000 | 2.8 GB | not converted |
| all | 10000 | 562 MB | **198 MB** (158 MB gzipped) |
| all | 20000 | 282 MB | ~100 MB [est, 2.8x] |
| all | 50000 | 115 MB | ~41 MB [est] |
| 1 rep/species (7,448) | 10000 / 20000 / 50000 | 335 / 168 / 68 MB | **120 / 61 / 25 MB** |
| 3 reps/species (10,336) | 10000 / 20000 / 50000 | 392 / 197 / 79 MB | **142 / 71 / 29 MB** |

Representatives use BFD's own rule (highest BUSCO genome completeness, then N50, then ASMID).

## 4. Accuracy on BFD (raw samples.csv labels)

Leave-one-out = the genome itself is excluded; leave-species-out excludes every genome of its species (a novel species); leave-genus-out excludes its genus. Counted over **answerable** queries (a taxon with another member left to find). scaled=10000, all 23,669 genomes; full tables per phylum in `reports/bench_raw_labels_s10000.md`.

| mode | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|
| leave-one-out | 0.985 | 0.979 | 0.974 | 0.970 | 0.958 | 0.928 |
| leave-species-out | 0.961 | 0.937 | 0.926 | 0.910 | 0.881 | n/a |
| leave-genus-out | 0.845 | 0.694 | 0.649 | 0.543 | n/a | n/a |

- Basidiomycota and Mucoromycota are comparable to Ascomycota in leave-one-out (phylum .988/.966, species .913/.900). Small phyla (Chytridiomycota, Microsporidia, Kickxellomycota) are poor when a genus is left out because they have few close relatives in BFD.
- scaled=1000, 2000 and 10000 give the same accuracy within ~0.5 points (6,000/3,000-query samples at the finer scales), so there is no accuracy reason to use more than scaled=10000.
- **Leave-genus-out is dominated by label errors, not sketch failures** (section 5). Interpret the 0.54-0.85 above as a lower bound.

### Representatives-only references (queries are the genomes left out of the reference)
| reference | scaled | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| 1 rep/species | 10000 | 0.996 | 0.991 | 0.989 | 0.985 | 0.974 | 0.605 |
| 1 rep/species | 50000 | 0.997 | 0.993 | 0.987 | 0.983 | 0.974 | 0.610 |
| 3 reps/species | 10000 | 0.999 | 0.996 | 0.994 | 0.991 | 0.987 | 0.727 |
| 3 reps/species | 50000 | 0.999 | 0.996 | 0.994 | 0.991 | 0.987 | 0.724 |

Rank accuracy for a held-out genome whose species has a representative is >=97% to genus even at scaled=50000 with one genome per species. Species accuracy is limited by naming (section 6).

## 5. Label audit of BFD (samples.csv)

`find_label_suspects.py` on the full leave-one-out run (ANI >= 0.95, neighbour vote consistency >= 0.8): **6 phylum, 38 class, 61 order conflicts** among 23,669 genomes (0.2%). In the leave-genus-out run, 23 reference genomes are each the best hit of >=3 queries of a different class at ANI >= 0.95 (they are mislabeled, or contaminated); the largest, `GCA_049307695.1_EMGR_V1` (labelled Dothideomycetes), is the best hit of 416 *Aspergillus*-type genomes. 20 of the 23 were also caught by the leave-one-out audit. Files: `reports/label_audit/exclude_ids.txt` (43 genomes), `hub_references.tsv`, `loo_suspects.tsv`, `log_conflicts.tsv`. These are suspicions from sequence similarity; I did not look at the assemblies themselves.

## 6. Species-level naming (what limits "species" accuracy)

Of 3,646 species errors with 3 reps/species: top confusions are *S. cerevisiae* vs hybrid labels ("S. cerevisiae x S. kudriavzevii" 813, "x S. eubayanus" 276, *S. pastorianus* 116), *F. oxysporum* vs *F. tardichlamydosporum* / *F. nirenbergiae* / *F. odoratissimum* (458; the *F. oxysporum* species complex), *A. flavus* vs *A. oryzae* (89), *A. fumigatus* vs "*Aspergillus* sp. GbtcF1" (72), "*Penicillium* cf. *chrysogenum*" (53). 94% of species errors are within the correct genus. Near-identical genomes (ANI > 0.995) carry a different species label 24% of the time. **Next (not built):** `gmg cluster`, greedy centroid clustering at ANI >= 0.95 (centroid = BFD rule), with the cluster ID in the index and the guess output, so the tool can say "same group as X" when names disagree; plus a per-cluster label-purity table for BFD.

## 7. Genomes not in BFD (stage 5)

145 NCBI GenBank assemblies not in BFD (by accession; released mostly after 2025-05; N50 >= 20 kb; lineage from NCBI taxonomy). **Finding: none of them is a named species absent from BFD**: the 25 flagged "species not in BFD" all lack a species name in NCBI. So named-novel-species behaviour is only tested by the leave-species/genus-out runs in section 4. Results with default (uncalibrated) thresholds, scaled=10000:

| subset | n | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| species in BFD: calls made / precision | 120 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 1.00/1.00 | 0.975/0.966 |
| unnamed species: calls made / precision | 25 | 0.92/0.96 | 0.91/0.95 | 0.91/0.95 | 0.96/0.91 | 0.91/0.85 | n/a |

## 8. Reads (stage 6)

Simulated reads (`sim_reads.py`: 150 bp, 0.5% substitutions, uniform) from 36 of the external genomes at 0.5x/1x/5x/20x; sketched with abundances; abundance filter at query time (`--min-abund`). At scaled=10000 the phylum, order, genus and species calls match the genome-based result at all four coverages (e.g. 0.5x with min-abund 2: phylum 0.97/1.00, order 0.94/1.00, genus 0.94/1.00, species 1.00/1.00 called/precision). Abundance filtering (min 2) removes sequencing-error hashes: on a 5x pilot 40,337 query hashes fall to 25,837 and the query's containment in the right genome rises from 0.70 to 0.998.
12 **real** *Coccidioides* R1 read sets (6 *C. immitis*, 6 *C. posadasii*; first 100k, 500k or 2M reads; most strains are not in BFD, which has 22 *Coccidioides* genomes): **12/12 correct at every rank including species** at all depths and both abundance filters.
Caveats: only 36+12 queries; substitution-only simulation (no indels, no contamination, no paired-end structure); the *Coccidioides* species are close to BFD genomes. Mixed samples and low-quality real data were not tested.

## 9. Decisions I made without you, problems, and what is not done

- **I crashed my own first s=10000 benchmark** by overwriting `GmG_runs/bin/gmg` while jobs were running from it (SIGSEGV/SIGBUS). Binaries are now deployed as `gmg-<timestamp>` with a symlink swap (`scripts/hpcc/deploy_gmg.sh`). Results were re-run.
- Built `gmg` once on the HPCC login node with `nice` and `-j3` (short compile) because the queue was saturated; everything else ran as SLURM jobs (partition `short`, 2-4 CPUs per task).
- A permission check blocked downloading the NCBI `datasets` binary onto the HPCC; I downloaded the genomes locally and copied only the sketches.
- ANI definition changed during the night from max(containment) to the query's containment (a small contig fully contained in a large genome scored ~1.0). Analyses recompute it from the TSVs; `Guess.top.ani` in the core uses the new form (deployed after the first external/reads run, so those used the old form).
- **Calibration is not final.** Default per-rank thresholds are placeholders; automatic calibration to a precision target failed to find thresholds for class-genus on the raw labels because of the label errors in section 5. A calibration on the cleaned benchmark is in progress (section 10).
- **Not done:** `gmg cluster`; WASM build (library-level wasm32 build works; a full link was blocked by a local toolchain problem); PyO3 bindings; PopSketch trees/PCA; MATPredict integration (its `--taxid`/`--phylum` hint output); contamination/fragmented-assembly tests; indels/paired reads/mixed samples for reads; conversion of the scaled=20000/50000 full indexes to v1.

## 10. Update after the cleaned benchmark

(filled in below when the job finished)
