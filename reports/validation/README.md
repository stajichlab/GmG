# Validation summary (2026-10-28)

What was tested, on data the index was not built from, and what it shows. Scripts: `scripts/validate_run.py` (scoring), `heldout_calibration.py`, `novelty_curve.py`.
Details per dataset are in the linked files. All numbers are from runs of `gmg run` with the BFD index (23,669 genomes, scaled = 10000) and the default `known` thresholds.

| # | dataset | n | truth | file |
|---|---|---|---|---|
| 1 | NCBI fungal assemblies released after the BFD freeze, not in BFD, **no quality filter** | 1,257 | NCBI taxonomy (only 52 have a species name; 1,240 a genus) | `post_freeze.md` |
| 2 | Real Illumina reads, *Coccidioides* (R1 only), first 100k / 500k / 2M reads | 558 each | corrected species of the MATPredict manifest | `cocci_reads_n*.md` |
| 3 | Held-out calibration: thresholds fitted on half the genera, tested on the other half, and fitted on Ascomycota, tested on other phyla | 70,878 bench queries | BFD labels (43 suspect genomes removed) | `heldout_calibration.md` |
| 4 | False-call rate for taxa absent from the reference (leave-genus-out / leave-species-out) | 23,6k each | BFD labels | `novelty_curve.md` |
| 5 | Real RNA-seq reads (BFD, species-labelled, non-empty files only) | 345 | BFD species lineage | running |

## What works

* **Independent assemblies.** 1,257 genomes, all processed (0 errors). Precision of the calls made: phylum 0.990, class 0.991, order 0.986, family 0.983, genus 0.971; calls made for 95% (phylum) to 90% (genus). Chromosome-level
  assemblies: phylum 1.000, genus 0.979 at 99% called; N50 ≥ 1 Mb: phylum 0.997, genus 0.979. The 52 genomes with a named species (all in BFD): species 1.000 precision, 75% called.
* **Real reads.** 558 *Coccidioides* read sets: genus called for ≥ 99.6% and correct for 99.1-99.3%, **even from 100k reads** (about 0.5x coverage, 160 sketch hashes). *C. immitis* versus *C. posadasii*: species called for
  88-93% with precision 0.996 at every depth; the two relabelled strains come out as their corrected species.
* **Calibration generalises to unseen genera.** With thresholds fitted on one half of the genera, the other half has precision 0.945-0.962 against 0.950-0.965 in-sample.

## What does not (yet)

1. **A taxon that is not in the reference still gets a name.** With the default genus cutoff (ANI ≥ 0.7765) a genus is named for **58%** of queries whose genus was left out of the reference (BFD, n = 23k), and 13 of 22 (**59%**) of the independent genomes
   from genera absent from BFD; their family, order, class and phylum calls are also weaker (precision 0.74, 0.73, 0.74, 0.80 for family to phylum, n = 39). The cutoff was fitted where relatives exist and never saw this case.
   Raising the genus cutoff to **ANI ≥ 0.90** lowers the false-call rate to 18% (BFD) / 14% (independent) and *raises* precision on known genera (0.985 vs 0.969) for 8 points of coverage (0.986 → 0.906) (`novelty_curve.md`;
   the independent set agrees but was not used to choose it). Species: at ANI ≥ 0.9938 a species is still named for 11.7% of queries whose species is absent. File: `data/thresholds/bfd_s10000.balanced.json` (not the default).
2. **Thresholds fitted on Ascomycota under-deliver elsewhere.** Tested on the other phyla precision falls to 0.87-0.93 (genus 0.873, species 0.892, class 0.927); on the independent set Basidiomycota genus precision is 0.942 [0.893, 0.969] and Microsporidia phylum 0.889.
   The family threshold is unstable between splits (coverage 0.40-0.98). The `novel` thresholds do not generalise beyond phylum and class.
3. **Tiny or fragmented inputs get no call, correctly, but are described wrongly.** 59 of the 1,257 assemblies have N50 < 10 kb; they are partial assemblies (median 1.5 Mb, 148 sketch hashes) and 48 are `no_match`
   (including 23 of 24 Neocallimastigomycota, median 2 Mb). QC flags (`short_genome`, `few_hashes`) say so, but the report text says "not a fungus, or a lineage the reference does not cover".
4. **Reads: more reads, fewer species calls.** From 100k to 2M reads the share of *Coccidioides* read sets that clear the species ANI cutoff falls from 92.1% to 88.5% (median ANI 0.9968 → 0.9957): repeated sequencing errors
   survive the abundance filter and lower the containment ANI. Reads need their own coverage-aware calibration (#32).
5. **Limits of the evidence.** Species-level truth is thin in the independent set (52 genomes); the novel-genus set is small (22 with a named genus); labels are NCBI taxonomy and BFD metadata (which has known errors, `reports/label_audit/`);
   the *Coccidioides* test is two closely related species of one genus and some strains may have genomes in BFD (22 *Coccidioides* genomes); calibration used the same BFD benchmark data for fitting and the held-out *splits*, not an external set.

## Follow-ups (issues)

default genus/species policy for absent taxa; clade-aware calibration; tiny-input messaging; reads calibration (#32); RNA-aware index (#33); error codes (#55); stress tests (#41); more independent sets (#40).
