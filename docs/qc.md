# Input QC flags

Flags go to the `qc_flags` column of `results.tsv`, to `qc` in `results.jsonl` and to the report. **They explain a result;
they never change a call.**

| flag | rule | why this threshold |
|---|---|---|
| `short_genome` | assemblies: the sketch-based size estimate (hashes x scaled) is below the 1st percentile of the **called phylum** in the reference (0.5th percentile of all reference genomes when no phylum is called or the phylum has fewer than 50 genomes) | the genome is probably partial |
| `long_genome` | same estimate above the 99th percentile of the phylum (99.5th of all) | contamination, a concatenated assembly, or a lineage the reference does not cover |
| `fragmented` | N50 below 10 kb | below this the assembly is draft-quality; the sketch is usually still informative, a flagged input that comes out `no_match` is more likely partial than non-fungal (see #59, `reports/validation/README.md`) |
| `high_N` | more than 5% of the bases are not A/C/G/T | gaps or masked sequence lower the number of hashes |
| `few_hashes` | fewer than 100 sketch hashes (at scaled 10000: under ~1 Mb of sequence) | a conservative floor: in the completeness stress test inputs with a median of 453 hashes (10% of a genome) still got correct genus calls (19/19 called), but nothing was tested below that, and at 100 hashes a single chance match is a large share of the evidence |
| `identical_to_reference` | assemblies: ANI >= 0.9999 and the reference is contained >= 0.999 | the input is (a copy of) a genome in the reference: the call is a lookup, not an inference |
| `not_fungal` | the phylum call or the best hit is a non-fungal decoy genome, with best-hit ANI >= 0.90 | the input has a close relative among the non-fungal decoys (plants, animals, protists): most likely not a fungus; any fungal rank is not meaningful |
| `distant_non_fungal_hit` | the nearest genomes are non-fungal decoys, but the best-hit ANI is below 0.90 | no close relative anywhere in the reference: a divergent fungus and a divergent non-fungus both look like this (best hit by chance, ANI 0.75-0.88). Not a verdict. The call is made again without the decoys, so a divergent fungus gets the call it would get from an index without them (never a Chordata or Streptophyta phylum); a divergent non-fungus gets a low-confidence fungal call or `no_match`, so treat the call as unknown |
| `weak_phylum` | assemblies: a phylum was called but the best-hit ANI is below 0.85 | the phylum is right only about 4 times in 5 for divergent genomes and is often a chance best hit (RefSeq truth set: phylum precision 0.998 for leave-one-out and 0.981 for leave-genus-out at best-hit ANI >= 0.85, 0.858 for leave-genus-out overall). The call is kept; filter on the flag if you need the stricter set. Not set for reads (errors lower their ANI) |
| `mixture_suspected` | a genus other than the best hit's explains at least 1% of the query hashes (and at least 20 of them) that no genome of the best hit's genus has, at least 30% of everything the best hit's genus does not explain, and at least half of its hashes sit in one of its genomes (one organism, not conserved sequence spread over many). Reference exclusions and decoys are respected. Not set for `not_fungal` inputs. `explained_fraction`, `second_taxon`, `second_explained_fraction` give the numbers | 0 of 2,191 pure or contaminated (bacteria, human) assemblies were flagged (RefSeq truth set with its own genome, species or genus hidden, 24 hosts); 2 of 64 real read libraries (probable *Dimargaris* + *Cokeromyces* co-cultures). Detection needs the second genus in the reference: with it, about 2% of the bases in assemblies and 10-20% of a 4x read set (minority coverage 0.4-0.8x, the abundance cutoff hides less); a minority whose genus is not in the reference, or of the same genus, is never flagged. The fractions are lower bounds of the share of sequence content (reads: not of reads). Numbers: `reports/validation/mixture_detector.md` (#31). The call is never changed |
| `high_error_reads` | reads: the k-mer survival estimated from the abundance histogram, `coverage_est` / (bases / size of the nearest reference genome), is <= 0.5 (error rate about 3.3% or more, as (1-e)^21), and the base coverage is >= 1x | sequencing errors remove k-mers: ONT R9.4 and PacBio CLR (5-15%) come out with a genus but rarely a species call. Measured on simulated long reads (6 genomes, 0.5-15% error, 1-30x): 0 of 60 read sets flagged at 0.5% and 1% error, 4 of 30 at 2%, 29 of 30 at 5% and 6 of 6 from 5x at every higher error; in 20 real runs 11 of 11 at >= 7% measured error and 1 of 9 at <= 2.9%. No flag does not prove low error (no estimate below 1x or without a histogram peak). `qc.error_est` and `qc.kmer_survival_est` in `results.jsonl`; the call is unchanged (reports/validation/long_reads.md, #44) |

The size rule is a range from the reference, not a fixed number, because genome sizes differ by lineage. BFD (23,669
genomes, estimated size = hashes x scaled):

| group | n | p1 | p5 | median | p95 | p99 |
|---|---|---|---|---|---|---|
| all reference genomes | 23,669 | 6.5 Mb | 10.0 | 30.9 | 51.9 | 71.8 |
| Ascomycota | 19,403 | 8.5 | 10.2 | 31.2 | 51.0 | 59.3 |
| Basidiomycota | 3,275 | 7.2 | 11.3 | 30.2 | 64.2 | 236.8 |
| Microsporidia | 163 | 0.8 | 1.6 | 4.2 | 15.8 | 31.4 |
| Glomeromycota | 60 | 26.6 | 43.8 | 96.8 | 296.1 | 376.4 |

A fixed 5 Mb / 150 Mb rule (what earlier versions used) would flag half of all Microsporidia as short and the rusts
(Basidiomycota, up to ~240 Mb) as long; with the phylum range they are normal. For an unclassified input the fallback
range is 3.6-98.7 Mb (p0.5-p99.5 of all genomes). The report card shows "Typical size" for the called lineage.

**Checks.** The completeness stress test (`reports/validation/stress_tests.md`: 10-75% of a genome) and the independent
assemblies with N50 < 10 kb are the test cases; `crates/gmg-cli/tests/run.rs` has a unit-level test of the size rule
(`short_genome` for a 20 kb piece of a 150-250 kb lineage, `long_genome` for a 900 kb input, none for the whole genome).

**Not done:** flags for mixtures (no reliable signal, `reports/validation/stress_tests.md`), a contig-count flag (N50 covers
fragmentation) and per-lineage N50 expectations.
