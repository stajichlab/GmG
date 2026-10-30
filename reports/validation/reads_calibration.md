# Reads: coverage estimate and abundance cutoff (#32)

Index `bfd-2026.10.4-s10000.gmgi`; `gmg` built from branch `reads-calibration` (the code of PR #109). Scripts: `scripts/hpcc/32_reads_calibration.sbatch` (cutoff comparison),
`scripts/hpcc/33_coverage_sim.sbatch` (simulated reads of known coverage), `scripts/reads_calibration_split.py`, `scripts/reads_calibration_summary.py` (the tables below).
Per-sample output stays on the HPCC (`GmG_runs/reads/calibration`).

## Summary and decisions

1. **`qc.coverage_est` (new, schema 7; `coverage_est` column in results.tsv).** The genome coverage of the sketched reads, from the abundance histogram: the peak after the error valley, or below about 3x a Poisson
   ratio of the numbers of hashes seen 3 and 2 times (definition in `docs/results-schema.md`, code in `crates/gmg-core/src/coverage.rs`). It is **informational, not used for any call**. Accuracy (section B):
   * clean simulated reads, 1x-30x: median ratio 0.81-0.97 (0.90 at 0.5% errors), 89-100% of read sets within 30% of the truth; below 0.3x it is noisy and often absent (too few hashes seen twice);
   * it is biased low by the sequencing error: genomic k-mers lose a factor (1 - e)^k, about 0.90 at 0.5% and 0.81 at 1% errors (the 0.81 for the 1% simulations);
   * real clean libraries (558 *Coccidioides* read sets, truth = bases / a nominal 28.5 Mb genome): median ratio 0.84 at 500k reads and 0.76 at 2M reads (within 30%: 83% and 58%);
   * **ZyGoLife libraries do not agree** (median ratio 0.32-8.9, within 30% for 0-43%). Their histograms are over-dispersed (a big singleton spike and hundreds of hashes seen more than 12 times,
     with only about a quarter of the assembly's k-mers seen twice at ~15x of bases in the worst examples I inspected), so either the libraries are not what bases / assembly size assumes (contamination, duplicates, high-copy sequence) or
     the estimator fails on them; I could not tell which without mapping the reads. The median abundance of the kept hashes (`median_kept`) is no better there and is biased high at low coverage
     (median 4.4 times the true coverage for *Coccidioides* at 100k reads, about 1x). Read the number as "k-mer coverage of the dominant genome", not as bases / genome size.
2. **`--min-abund auto` (new, opt-in; `gmg run` and `gmg guess`): the cutoff at the valley of the histogram.** **Default stays 2.** Decision on the evaluation half (section A): auto is identical to 2 at low depth
   (the histogram has no valley: it chose 2 for 345 of 345 RNA sets and 558 of 558 *Coccidioides* sets at 100k reads) and at depth (about 3-10x) it picked 3-20 for roughly half of the sets, with changes of a few
   calls in both directions on the evaluation half (RNA species at 2M reads +4 calls, 0.68 -> 0.70; *Coccidioides* genus and species at 500k/2M reads -1 to -2 calls; ZyGoLife genus at 1.5M reads -1 call
   but phylum precision 0.918 -> 0.938). The rule for making it the default was "no loss of genus or species calls or precision"; it loses 1-2 calls in 6 of the 12 evaluation cells at the two deepest depths
   of each collection (4 *Coccidioides*, 2 ZyGoLife; none in RNA-seq), so it is not the default.
   **Fixed 1 is unusable for reads** (species called for 8-10% of *Coccidioides* sets at 2M reads instead of 96%; ZyGoLife phylum precision 0.81-0.83 instead of 0.88-0.89 at 300k-1.5M reads on the calibration half):
   errors are kept. The fixed cutoff 2 stays.
3. **Paired reads.** R1/R2 and lane files of one sample are merged before sketching, with abundances summed. Test `paired_and_lane_split_reads_equal_the_concatenated_file` now checks the equality with the default
   cutoff and with `auto` too, for status, call, `n_hashes`, `shared`, `total_bases` and `coverage_est`. Paired against single-end on real reads: a merged pair is worth about twice the reads of R1 alone
   (`reads_titration.md`).
4. **Calibration and evaluation sets.** Every sample id goes to the calibration half (`cal`: first byte of SHA-1(id) even) or the evaluation half (`eval`: odd), the same way whatever else is in the table
   (`scripts/reads_calibration_split.py`). The choice between fixed 1, fixed 2 and auto was made on the calibration half and checked on the evaluation half; both halves tell the same story. The cutoff has no fitted
   parameter other than the valley rule, so the split guards against choosing by the numbers that are then reported, not against over-fitting.

## A. Abundance cutoff: fixed 1, fixed 2, auto

Cells: share of inputs with a call / precision of the calls (calls/inputs). `cal` = calibration half, `eval` = evaluation half (SHA-1 of the sample id).


**ZyGoLife DNA reads (genus truth = lab label)**

| reads per file | cutoff | cal phylum | cal genus | eval phylum | eval genus |
|---|---|---|---|---|---|
| 15000 | 1 | 0.79 / 0.873 (55/70) | 0.73 / 0.745 (51/70) | 0.73 / 0.944 (36/49) | 0.67 / 0.758 (33/49) |
| 15000 | 2 | 0.76 / 0.887 (53/70) | 0.76 / 0.736 (53/70) | 0.73 / 0.944 (36/49) | 0.69 / 0.735 (34/49) |
| 15000 | auto | 0.76 / 0.887 (53/70) | 0.76 / 0.736 (53/70) | 0.73 / 0.944 (36/49) | 0.69 / 0.735 (34/49) |
| 75000 | 1 | 0.87 / 0.885 (61/70) | 0.77 / 0.704 (54/70) | 0.94 / 0.891 (46/49) | 0.71 / 0.771 (35/49) |
| 75000 | 2 | 0.86 / 0.883 (60/70) | 0.79 / 0.709 (55/70) | 0.88 / 0.954 (43/49) | 0.76 / 0.757 (37/49) |
| 75000 | auto | 0.86 / 0.883 (60/70) | 0.79 / 0.709 (55/70) | 0.88 / 0.954 (43/49) | 0.76 / 0.757 (37/49) |
| 300000 | 1 | 0.97 / 0.809 (68/70) | 0.74 / 0.750 (52/70) | 1.00 / 0.878 (49/49) | 0.71 / 0.771 (35/49) |
| 300000 | 2 | 0.91 / 0.891 (64/70) | 0.77 / 0.741 (54/70) | 0.98 / 0.938 (48/49) | 0.78 / 0.763 (38/49) |
| 300000 | auto | 0.91 / 0.891 (64/70) | 0.77 / 0.741 (54/70) | 0.98 / 0.958 (48/49) | 0.78 / 0.763 (38/49) |
| 1500000 | 1 | 0.99 / 0.826 (69/70) | 0.71 / 0.760 (50/70) | 0.98 / 0.917 (48/49) | 0.71 / 0.771 (35/49) |
| 1500000 | 2 | 0.99 / 0.884 (69/70) | 0.79 / 0.746 (55/70) | 1.00 / 0.918 (49/49) | 0.76 / 0.784 (37/49) |
| 1500000 | auto | 0.99 / 0.913 (69/70) | 0.79 / 0.746 (55/70) | 0.98 / 0.938 (48/49) | 0.73 / 0.778 (36/49) |

Cutoffs chosen by `auto` (cutoff: number of read sets): 15000 reads: 2: 120; 75000 reads: 2: 120; 300000 reads: 2: 112, 4: 1, 7: 2, 8: 1, 9: 1, 10: 1, 13: 1, 18: 1; 1500000 reads: 2: 68, 3: 14, 4: 7, 5: 5, 6: 4, 7: 4, 8: 5, 9: 2, 10: 5, 13: 1, 14: 1, 18: 1, 20: 2, 24: 1


**Coccidioides DNA reads (species truth)**

| reads per file | cutoff | cal genus | cal species | eval genus | eval species |
|---|---|---|---|---|---|
| 100000 | 1 | 1.00 / 0.992 (258/258) | 0.66 / 1.000 (170/258) | 1.00 / 0.997 (299/300) | 0.64 / 1.000 (192/300) |
| 100000 | 2 | 1.00 / 0.992 (258/258) | 0.97 / 0.992 (249/258) | 0.99 / 0.997 (298/300) | 0.94 / 0.996 (281/300) |
| 100000 | auto | 1.00 / 0.992 (258/258) | 0.97 / 0.992 (249/258) | 0.99 / 0.997 (298/300) | 0.94 / 0.996 (281/300) |
| 500000 | 1 | 1.00 / 0.992 (258/258) | 0.39 / 1.000 (100/258) | 1.00 / 0.997 (299/300) | 0.41 / 1.000 (123/300) |
| 500000 | 2 | 1.00 / 0.992 (258/258) | 0.96 / 0.996 (247/258) | 1.00 / 0.993 (300/300) | 0.96 / 0.997 (287/300) |
| 500000 | auto | 1.00 / 0.992 (257/258) | 0.96 / 0.996 (247/258) | 0.99 / 0.993 (298/300) | 0.95 / 0.997 (286/300) |
| 2000000 | 1 | 1.00 / 0.992 (258/258) | 0.10 / 1.000 (26/258) | 1.00 / 0.993 (300/300) | 0.08 / 1.000 (24/300) |
| 2000000 | 2 | 1.00 / 0.992 (258/258) | 0.96 / 0.996 (248/258) | 1.00 / 0.997 (299/300) | 0.96 / 0.997 (289/300) |
| 2000000 | auto | 1.00 / 0.992 (258/258) | 0.97 / 0.996 (250/258) | 0.99 / 0.997 (297/300) | 0.96 / 0.997 (288/300) |

Cutoffs chosen by `auto` (cutoff: number of read sets): 100000 reads: 2: 558; 500000 reads: 2: 555, 7: 1, 9: 1, 12: 1; 2000000 reads: 2: 210, 3: 250, 4: 78, 5: 16, 6: 1, 9: 1, 16: 1, 17: 1


**BFD RNA-seq (species truth)**

| reads per file | cutoff | cal genus | cal species | eval genus | eval species |
|---|---|---|---|---|---|
| 100000 | 1 | 0.92 / 0.994 (167/182) | 0.60 / 0.955 (111/184) | 0.96 / 1.000 (153/160) | 0.70 / 0.964 (112/161) |
| 100000 | 2 | 0.91 / 0.994 (166/182) | 0.64 / 0.957 (117/184) | 0.95 / 1.000 (152/160) | 0.73 / 0.966 (118/161) |
| 100000 | auto | 0.91 / 0.994 (166/182) | 0.64 / 0.957 (117/184) | 0.95 / 1.000 (152/160) | 0.73 / 0.966 (118/161) |
| 500000 | 1 | 0.92 / 0.994 (167/182) | 0.55 / 0.960 (101/184) | 0.96 / 1.000 (153/160) | 0.63 / 0.971 (102/161) |
| 500000 | 2 | 0.91 / 0.994 (166/182) | 0.65 / 0.958 (119/184) | 0.96 / 1.000 (153/160) | 0.73 / 0.966 (117/161) |
| 500000 | auto | 0.91 / 0.994 (166/182) | 0.65 / 0.958 (119/184) | 0.96 / 1.000 (153/160) | 0.73 / 0.966 (117/161) |
| 2000000 | 1 | 0.92 / 0.994 (168/182) | 0.40 / 0.959 (73/184) | 0.96 / 1.000 (153/160) | 0.43 / 0.986 (70/161) |
| 2000000 | 2 | 0.93 / 0.994 (169/182) | 0.59 / 0.963 (109/184) | 0.95 / 1.000 (152/160) | 0.68 / 0.963 (109/161) |
| 2000000 | auto | 0.93 / 0.994 (170/182) | 0.61 / 0.955 (112/184) | 0.95 / 1.000 (152/160) | 0.70 / 0.965 (113/161) |

Cutoffs chosen by `auto` (cutoff: number of read sets): 100000 reads: 2: 345; 500000 reads: 2: 342, 3: 2, 10: 1; 2000000 reads: 2: 174, 3: 21, 4: 46, 5: 46, 6: 25, 7: 21, 8: 4, 9: 2, 10: 1, 11: 4, 13: 1


## B. Accuracy of `coverage_est`

Ratio = estimated / true coverage; `within 30%` = share of read sets with ratio in 0.7-1.3.


**Simulated reads (production scaled 10000), default cutoff 2**

| set | n | median ratio | 10-90% range | within 30% |
|---|---|---|---|---|
| 0.1x (3 genomes x 3 error rates) | 3 | 0.70 | 0.60-0.90 | 0.67 | no estimate: 6
| 0.3x (3 genomes x 3 error rates) | 9 | 0.80 | 0.57-2.00 | 0.33 |
| 1x (3 genomes x 3 error rates) | 9 | 0.97 | 0.86-1.24 | 1.00 |
| 3x (3 genomes x 3 error rates) | 9 | 0.93 | 0.77-0.99 | 1.00 |
| 10x (3 genomes x 3 error rates) | 9 | 0.81 | 0.69-1.04 | 0.89 |
| 30x (3 genomes x 3 error rates) | 9 | 0.92 | 0.77-1.00 | 1.00 |
| error 0.1% (all coverages) | 16 | 0.98 | 0.90-1.53 | 0.88 |
| error 0.5% (all coverages) | 16 | 0.90 | 0.77-1.22 | 0.94 |
| error 1% (all coverages) | 16 | 0.81 | 0.60-0.93 | 0.69 |

**Real reads, default cutoff 2: coverage = sketched bases / genome size**

| set | n | median ratio | 10-90% range | within 30% |
|---|---|---|---|---|
| zygo, first 15000 reads per file (true coverage from its assembly) | 30 | 8.93 | 3.70-17.75 | 0.00 |
| zygo, first 75000 reads per file (true coverage from its assembly) | 43 | 2.67 | 0.60-5.44 | 0.09 |
| zygo, first 300000 reads per file (true coverage from its assembly) | 46 | 0.73 | 0.50-1.56 | 0.43 |
| zygo, first 1500000 reads per file (true coverage from its assembly) | 46 | 0.32 | 0.13-2.27 | 0.11 |
| cocci, first 100000 reads per file (true coverage from a nominal 28.5 Mb genome) | 553 | 1.11 | 0.61-2.04 | 0.49 |
| cocci, first 500000 reads per file (true coverage from a nominal 28.5 Mb genome) | 558 | 0.84 | 0.65-0.99 | 0.83 |
| cocci, first 2000000 reads per file (true coverage from a nominal 28.5 Mb genome) | 558 | 0.76 | 0.46-0.88 | 0.58 |

## Not covered

* The error bias of the estimate is documented, not corrected (an error-rate estimate from the singleton spike would allow it).
* No mapping of the ZyGoLife reads to their assemblies, so the disagreement of section B on those sets is explained, not proven.
* `auto` was compared with 1 and 2 on three collections at 3 to 4 depths each; sets with a deep, clean library (above 20x) are only represented by the simulations. If a deep clean collection shows a gain
  from the valley rule, `--min-abund auto` is the switch to use (`gmg run` and `gmg guess`), and the default can be revisited with that evidence.
* Single-end reads below 100 bp and long reads (#44).
