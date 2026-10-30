# Validation summary (2026-10-28)

What was tested, on data the index was not built from, and what it shows. Scripts: `scripts/validate_run.py` (scoring), `heldout_calibration.py`, `novelty_curve.py`.
Details per dataset are in the linked files. All numbers are from runs of `gmg run` with the BFD index (23,669 genomes, scaled = 10000) and the default `known` (balanced: genus ANI ≥ 0.90, tentative genus ≥ 0.7765) thresholds; the earlier permissive results are kept as `*_permissive.*`.

| # | dataset | n | truth | file |
|---|---|---|---|---|
| 1 | NCBI fungal assemblies released after the BFD freeze, not in BFD, **no quality filter** | 1,257 | NCBI taxonomy (only 52 have a species name; 1,240 a genus) | `post_freeze_balanced.md` |
| 2 | Real Illumina reads, *Coccidioides* (R1 only), first 100k / 500k / 2M reads | 558 each | corrected species of the MATPredict manifest | `cocci_reads_n*_balanced.md` |
| 3 | Held-out calibration: thresholds fitted on half the genera, tested on the other half, and fitted on Ascomycota, tested on other phyla | 70,878 bench queries | BFD labels (43 suspect genomes removed) | `heldout_calibration.md` |
| 4 | False-call rate for taxa absent from the reference (leave-genus-out / leave-species-out) | 23,6k each | BFD labels | `novelty_curve.md` |
| 5 | Real RNA-seq reads (BFD, species-labelled, non-empty files only) | 345 | BFD species lineage | `rna_reads_n*.md` |
| 6 | Independent lab and public assemblies: IndigoAg (199), ZyGoLife aaftf/sort (552), KassonLab (36), Antarctic MAGs (266, classification only), CCFEE AAFTF vs shovill (67 strains); registry `data/validation/datasets.tsv` | 1,120 | provider names / NCBI taxids / agreement | `external_summary.md` |

## Before and after the balanced default (#57, PR #60)

| | permissive (genus ≥ 0.7765) | balanced (genus ≥ 0.90) |
|---|---|---|
| Independent assemblies, genus: called / precision | 0.900 / 0.971 | 0.839 / **0.988** |
| Genus absent from BFD (n = 22): genus named | 59% (13), all wrong | **14% (3)**, all wrong |
| Inputs that get a `tentative_genus` instead | n/a | 78 of 1,257 |
| *Coccidioides* reads, genus called / precision (100k and 2M) | .996-.998 / .995 | unchanged |
| *Coccidioides* reads, species called / precision (2M) | .884 / .996 | unchanged |
| RNA-seq, genus called / precision at 100k reads | .956 / .997 | .930 / .997 |
| RNA-seq, genus called / precision at 2M reads | .985 / .994 | .939 / .997 |
| RNA-seq, species called / precision at 2M reads | .542 / .968 | .542 / .968 |

The change does what the novelty analysis predicted: false genus names for absent genera drop from 59% to 14%, genus precision rises to 0.988, and the cost is 6 points of genus coverage on the independent set (more on divergent RNA-seq: 4-5 points). Reads of
*Coccidioides* and all species calls are unaffected. The weak-genus inputs now carry a `tentative_genus` and status `low_confidence`, which is never a call.

## What works

* **Independent assemblies.** 1,257 genomes, 0 errors. Precision of the calls made: phylum 0.990, class 0.991, order 0.986, family 0.983, genus 0.988. The 52 genomes with a named species (all in BFD): species precision 1.000 at 75% called.
* **Real reads.** 558 *Coccidioides* read sets: genus called for ≥ 99.6% and correct for 99.1-99.3% even from 100k reads (about 0.5x coverage); species called for 88-93% with precision 0.996 at every depth. 345 RNA-seq sets: genus precision 0.997 at every depth.
* **Calibration generalises to unseen genera.** Thresholds fitted on one half of the genera give precision 0.945-0.962 on the other half against 0.950-0.965 in-sample.

## What does not (yet)

1. **Absent taxa still sometimes get a name.** 3 of 22 independent absent-genus genomes are still named (BFD leave-genus-out: 18%); species at ANI ≥ 0.9938 is named for 11.7% of queries whose species is absent. The 22-genome set is too small to confirm the rate (CI 0-0.56 for precision).
2. **Thresholds fitted on Ascomycota under-deliver elsewhere.** Held-out other phyla precision 0.87-0.93; Basidiomycota genus 0.942 on the independent set. Clade-aware cutoffs: #58.
3. **Tiny or fragmented inputs** (59 assemblies with N50 < 10 kb; 48 `no_match`, including 23 of 24 Neocallimastigomycota) get no call, correctly, but the report text blames the lineage rather than the input (#59).
4. **Reads: errors lower the species ANI.** Fixed in part (#32): read sets now use a species cutoff of 0.99 (assemblies 0.9938). *Coccidioides* species called goes from 91.4 / 92.8 / 88.4% to 95.0 / 95.7 / 96.2% (100k / 500k / 2M reads) at unchanged precision (0.994-0.996), and no longer falls with depth; RNA-seq species called 61 / 59 / 54% → 68 / 68 / 63% (precision 0.966-0.968, was 0.968-0.972). The cutoff was checked with *C. immitis* removed from the reference: its reads match *C. posadasii* at ANI 0.976-0.980 and are falsely named at 0.99 in 0.6% (1/175) of sets (10% at 0.98). One species pair only; RNA-seq still needs a transcript-aware reference (#33). Files: `*_reads099.*`, `cocci_absent_species.md`.
5. **Limits of the evidence.** Species truth is thin (52 genomes); labels are NCBI taxonomy and BFD metadata with known errors (`reports/label_audit/`); *Coccidioides* is two close species of one genus; calibration used BFD data for both fitting and the held-out splits.

## Follow-ups (issues)

clade-aware calibration (#58); tiny-input messaging (#59); reads calibration (#32); RNA-aware index (#33); error codes (#55); stress tests including contamination (#41); more independent sets (#40).
