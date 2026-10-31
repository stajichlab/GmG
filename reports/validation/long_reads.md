# Long reads (ONT, PacBio) and the error rate (#44)

Does 21-mer containment survive 5-15% error? Two parts: a **simulated sweep** (24 chromosome-level genomes x 8 error rates x 7 depths) and **20 public ONT / PacBio runs** whose genome is in the reference,
both classified with the own genome (or its species, or its genus) hidden from the reference (#85). Scripts: `scripts/sim_long_reads.py`, `long_reads_plan.py`, `long_reads_measure.py`, `long_reads_harvest.py`,
`long_reads_score.py`, `long_reads_tables.py`, `long_reads_flag_summary.py`, `scripts/hpcc/35_..39_long_reads_*.sbatch`; lists: `data/validation/long_reads_{sim_genomes,candidates,pilot,runs}.tsv`; results:
`long_reads_sim.tsv` (accuracy per tag x error x depth x rank), `long_reads_real.tsv` (one row per run and setting), `long_reads_flag_check.tsv`.

## Summary and decision

* **Errors withhold calls; they do not make wrong ones.** Genus precision is 1.00 at every error rate (0.5-15%) and depth in the sweep (own genome hidden), phylum precision >= 0.95, and with the genus hidden the wrong genus calls
  stay at 0-12% (4-9% typical) at every error rate, as for assemblies and short reads. In the 20 real runs every genus call was right (20/20, error 0.1-13%), phylum 20/20.
* **What errors cost is depth.** With the default (min-abund 2) the genus call reaches its plateau (0.75-0.79 of the genomes; the rest are withheld, mostly `ani_below_cutoff`: the nearest remaining genome is below the genus cutoff ANI 0.90) from 0.5x at <= 2% error, from 1x at 5-8%, from 5x at 10-15%.
  Species (cutoff ANI 0.99 for reads) is called for 0.82 of the genomes that have a conspecific in the reference at <= 5% error and 1-5x (0.71-0.82), and drops to <= 0.5 above 8% error at any depth. Real runs agree: species right for 6 of 20 (2 of 5 HiFi, 4 of 5 R10-like, 0 of 10 at >= 7% error), genus 20/20.
* **The "10% reads look like a 0.90 relative" worry holds only for `--min-abund 1`.** With min-abund 1 the median ANI is 0.896 at 10% error (0.846 at 15%) and falls further with depth (error k-mers accumulate: 0.846 at 10%/30x), the genus is lost at >= 12% error (0.46 or less at 10%, 0 from 5x)
  and the species at >= 1% error (species recall 0.06 at best; at 0.5% error from 2x on). With the default min-abund 2 the error k-mers (nearly all seen once) are dropped, the ANI of the k-mers that remain stays 0.97-0.99, and the cost is that true k-mers
  also drop out: a k-mer survives with probability (1-e)^21 (0.11 at 10%) and is then seen on average depth x (1-e)^21 times.
* **min-abund 2 breaks at high depth, `auto` repairs it.** At 30x the error k-mers recur often enough to pass the filter: species recall 0.65 (0.5% error), 0.29 (1%), 0.00 (>= 2%), ANI 0.99 -> 0.96. `--min-abund auto` (valley of the histogram) gives 0.82 / 0.82 / 0.82 / 0.76 at 30x for 0.5 / 1 / 2 / 5% error and is
  identical to 2 elsewhere (it picks 2 when the histogram has no peak); above 8% error nothing recovers the species. No real run reached that depth range (200 Mb of a 9-40 Mb genome is 5-22x, 74x for the 2.7 Mb *E. hellem*), so real runs show auto = 2.
* **The error rate can be read from `coverage_est`.** `coverage_est` is the k-mer coverage of the error-free k-mers, so `coverage_est / (bases / genome size)` = (1-e)^21: 0.9 at 0.5% error, 0.8 at 1%, 0.65 at 2%, 0.35 at 5%, 0.12 at 10%
  (median over the sweep, tables below); the genome size is available as the size of the nearest reference genome (hashes x scaled). Short reads (the #32 simulations, 150 bp, 0.1-1% error, 0.1-30x) give 0.57-1.0.
* **Decision: one small change, no second index.**
  1. **New flag `high_error_reads`** (reads; `qc_flags`, plus `qc.error_est` and `qc.kmer_survival_est` in `results.jsonl`): the survival ratio above is <= 0.5 (error about 3.3% or more) and the base coverage is >= 1x. Evidence below. It explains a result, never changes a call.
  2. **No error-aware species cutoff.** The species call is already conservative (precision 1.00 at every error rate); lowering the cutoff with the error would trade precision for a call the data barely support at >= 8% error.
  3. **No change of default for the abundance filter.** `--min-abund auto` is recommended for long reads (and any read set) at >= 10x; in the sweep it equals the default in all but 10 of 112 genus and species cells, is better in 9 (up to +0.82 species recall at 30x) and lower by one genome in one (0.75 against 0.79 genus, 0.5% error, 10x). Not made the default: the only measured gain is at depths no real run here reached, and #32 found no gain for short reads.
  4. **No smaller k / finer scaled, no second index.** The genus, which is what long reads at 8-15% error can give, is called with precision 1.00 from 5x at 15% error with the existing index; the species is the rank that needs a lower error, not a different index.
  Recommended use: docs/limitations.md and docs/tutorial.md.

## Simulation (part 1)

**Genomes** (`data/validation/long_reads_sim_genomes.tsv`): 24 chromosome-level or complete assemblies in the reference (NCBI `assembly_level`), 2.6-55 Mb, 5 phyla (Ascomycota 12, Basidiomycota 8, Microsporidia 2, Mucoromycota 1, Mortierellomycota 1), 22 classes; 17 have another genome of the species in the reference
(species recall is counted over those 17), all have a congener. **Reads:** `scripts/sim_long_reads.py` (numpy only; PBSIM3 and Badread are not installed on the HPCC): log-normal length (median 12 kb, 95% in 4.5-32 kb), random strand and start; per base an error with probability e, a substitution / insertion / deletion
in the ratio 0.4 / 0.3 / 0.3; e = 0.5, 1, 2, 5, 8, 10, 12, 15%; one 30x stream per (genome, error) cut into nested 0.25, 0.5, 1, 2, 5, 10, 30x files (6,720 read sets). The realised error was checked by minimap2 alignment (`long_reads_measure.py`: 1.0 / 9.4 / 13.6% measured for 1 / 10 / 15% nominal on a test genome;
the alignment counts a substitution next to an indel once, so it reads lower). Errors are uniform and independent along the read; real errors cluster in homopolymers and low-quality stretches, so real k-mer survival is somewhat higher than (1-e)^21. **Classification:** `gmg run --kind reads` with the full index of release 2026.10.4 (23,669 genomes, s10000; the reps index has one genome per species, so hiding the genome would hide the species), `known` thresholds,
own genome hidden (`--exclude-level self`) with min-abund 2 (default), 1 and auto; species and genus hidden with min-abund 2. Binary `gmg-readscal` (`--min-abund auto`, `coverage_est`; the verdicts do not depend on the later flag changes). "Called" = share of the 24 inputs with a call at the rank; precision over the calls made; recall over the attainable inputs.


## Real runs (part 2)

**Runs** (`data/validation/long_reads_runs.tsv`, listed before the download): from the ENA portal API, fungal (tax_tree 4751) WGS runs with `instrument_platform` OXFORD_NANOPORE or PACBIO_SMRT, `library_source` GENOMIC, single-ended, >= 250 Mb of bases, whose BioSample is the BioSample of a chromosome-level or complete assembly in the reference,
with the species named the same way by the ENA record, the NCBI organism name and the reference label (644 eligible runs). The chemistry is only a guess from the metadata (year, model, title), so 80 candidates (16 per guessed class: HiFi, CLR, ONT R10, ONT 2022-23, ONT <= 2021; `long_reads_candidates.tsv`) were **pilot-measured**: the first 4 MB of each
file (0.3 GB in total) aligned to the run's own assembly with minimap2 (`long_reads_pilot.tsv`; 80 of 80 gave an estimate). The 20 runs were then taken from the first 75 pilot results (5 per measured class, <= 1 per species): 20 species, 20 genera, Ascomycota 17, Basidiomycota 2 (*Cryptococcus deneoformans* and *Rhizoctonia solani*), Microsporidia 1 (*Encephalitozoon hellem*).
The guessed chemistry was often wrong (guess -> class measured on the 200 Mb: 2 of 5 guessed HiFi are 0.9-1.2%, 1 of 5 guessed CLR is 2.9% (the others 11-13%), ONT guessed from 2022-23 or 2019-21 spans 7-12%, and the two ONT runs guessed R10.4 (2025) measure 1.7% and 2.7%, not 1%), which is why the measured error, not the metadata, defines the classes below. No ONT run at about 1% was in the list; the 0.5-3% class is 2 ONT and 3 PacBio runs.
**Download:** the first 200 Mb of bases of each run, streamed over HTTPS (`long_reads_harvest.py`, 3 connections; the connection is closed at 200 Mb; Range-resume on a dropped connection: 3 resumes in 2 runs): **20 of 20 runs reached 200 Mb, 3.05 GB of compressed bytes read in total** (61-214 MB per run), against a cap of 10 GB; the `fastq-dump -X` fallback was not needed. The raw reads are kept in `long_reads/real/raw` on the HPCC.
**Truth** = the lineage of the chromosome-level genome of the same BioSample; that genome is hidden (`ref` column + `--exclude-level self`), or its species or genus. Measured error = alignment of 500 reads of the 200 Mb (a lower bound: unaligned reads and clipped ends are not counted). Binary: private build of the branch (`gmg` 36987f8 with `high_error_reads`), index and settings as in part 1.

| run | platform, year | species | measured error | coverage (200 Mb / genome) | genus | species | ANI | coverage_est | flag | min-abund 1: genus / species | auto: species | species hidden: genus | genus hidden: genus |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SRR32051149 | Sequel II 2025 | *Alternaria brassicae* | 0.1% | 5.9x | right | - | 0.9946 | 5.01 |  | right / - | - | right | - |
| SRR25662223 | Sequel II 2024 | *Gaeumannomyces avenae* | 0.2% | 4.8x | right | right | 0.9984 | 3.76 |  | right / - | right | right | - |
| SRR28482132 | Revio 2024 | *Colletotrichum fructicola* | 0.3% | 3.5x | right | - | 0.9975 | 3.39 |  | right / - | - | right | - |
| SRR21003521 | Sequel II 2022 | *Fusarium verticillioides* | 0.3% | 4.8x | right | - | 0.9994 | 3.85 |  | right / - | - | right | - |
| SRR24086010 | Sequel II 2025 | *Lodderomyces parapsilosis* | 0.3% | 15.5x | right | right | 0.9984 | 13.02 |  | right / - | right | - | - |
| SRR32357977 | Revio 2026 | *Candida albicans* | 0.9% | 13.7x | right | right | 0.9990 | 8.01 |  | right / right | right | right | - |
| SRR35605210 | Revio 2026 | *Arcopilus aureus* | 1.2% | 5.9x | right | right | 0.9928 | 3.42 |  | right / - | right | - | - |
| SRR35627297 | MinION 2025 | *Cryptococcus deneoformans* | 1.7% | 10.4x | right | right | 0.9979 | 7.03 |  | right / - | right | right | - |
| SRR32027809 | GridION 2025 | *Histoplasma ohiense* | 2.7% | 5.1x | right | right | 0.9974 | 3.19 | yes | right / - | right | right | - |
| SRR25033382 | Sequel 2023 | *Maudiozyma bulderi* | 2.9% | 13.8x | right | - | 0.9864 | 7.02 |  | right / - | - | right | - |
| SRR17853474 | MinION 2023 | *Encephalitozoon hellem* | 7.3% | 73.9x | right | - | 0.9169 | 29.24 | yes | - / - | - | - | - |
| DRR261529 | MinION 2021 | *Aspergillus luchuensis* | 7.5% | 5.6x | right | - | 0.9926 | 1.59 | yes | right / - | - | right | WRONG |
| ERR10089557 | MinION 2022 | *Torulaspora quercuum* | 9.2% | 19.4x | right | - | 0.9692 | 1.25 | yes | - / - | - | - | - |
| SRR12403034 | MinION 2021 | *Ogataea polymorpha* | 9.3% | 22.3x | right | - | 0.9663 | 1.03 | yes | - / - | - | right | - |
| SRR19837494 | Sequel 2023 | *Debaryomyces hansenii* | 11.0% | 15.4x | right | - | 0.9617 | 1.05 | yes | - / - | - | right | - |
| SRR23085835 | Sequel II 2024 | *Candidozyma vulturna* | 11.7% | 15.6x | right | - | 0.9647 | 0.99 | yes | - / - | - | right | - |
| SRR11278194 | Sequel 2021 | *Monilinia vaccinii-corymbosi* | 11.8% | 6.7x | right | - | 0.9031 | 0.64 | yes | - / - | - | - | - |
| SRR8987927 | MinION 2019 | *Pichia kudriavzevii* | 12.1% | 18.2x | right | - | 0.9668 | 1.21 | yes | - / - | - | - | - |
| SRR19543541 | PromethION 2022 | *Rhizoctonia solani* | 12.1% | 4.8x | right | - | 0.9806 | 1.36 | yes | right / - | - | - | - |
| SRR10197407 | PacBio RS II 2020 | *Penicillium polonicum* | 13.2% | 6.1x | right | - | 0.9820 | 0.57 | yes | - / - | - | right | - |


Per run, own genome hidden, min-abund 2 unless stated; "-" = not called. (`*Encephalitozoon hellem*` is a 2.7 Mb genome, hence 74x.)

| measured class | runs | genus right | species right | genus or species wrong | `high_error_reads` | min-abund 1: genus / species right | genus hidden: wrong genus calls | species hidden: genus right (wrong) |
|---|---|---|---|---|---|---|---|---|
| HiFi (< 0.5%) | 5 | 5 | 2 | 0 | 0 | 5 / 0 | 0 of 5 | 4 (0) |
| R10-like (0.5-3%) | 5 | 5 | 4 | 0 | 1 (2.7%) | 5 / 1 | 0 of 5 | 4 (0) |
| R9-like (3-11%) | 4 | 4 | 0 | 0 | 4 | 1 / 0 | 1 of 4 | 2 (0) |
| CLR-like (>= 11%) | 6 | 6 | 0 | 0 | 6 | 1 / 0 | 0 of 6 | 3 (0) |

**Real runs against the simulation.** (1) Precision: no wrong genus or species call in 20 runs (sim: 1.00); with the genus hidden 1 of 20 gets a (wrong) genus call, 5% (sim 4-12%); with the species hidden 13 of 20 get the right genus and none a wrong one. (2) Species recall falls from 4 of 5 at 0.9-2.9% error to 0 of 10 at >= 7%, as in the sweep (0.82 at <= 5%, <= 0.5 above 8%);
the 3 HiFi runs and 1 R9 run without a species call have ANI >= 0.9926 and are withheld by `low_agreement` of the near-best hits (conflicting labels in a dense clade, e.g. *Fusarium verticillioides* ANI 0.9994), which is not an error effect. (3) min-abund 1: genus right for 12 of 20 (all 20 with min-abund 2), species for 1 of 20 (6 of 20): the sweep's collapse is real
(median ANI 0.89 at 10% error). (4) `--min-abund auto` equals 2 on every run (the depths are 3.5-22x for 19 runs; the sweep's difference starts at >= 10x with <= 5% error, and 74x for *E. hellem* has 7% error, where nothing recovers the species). (5) ANI of the (own-genome-hidden) top hit: 0.986-0.999 for the 0.1-2.9% runs, 0.903-0.993 at 7-13% (the genus is called from ANI 0.90 up).
(6) Depth: the genus is right from 3.5x at every error rate here; the sweep's 5x-at-15% need is met (4.8-6.7x for the three runs at 11.8-13.2% with < 7x).

**The flag on real runs.** `high_error_reads` fires for all 11 runs with measured error >= 7.3% and for 1 of the 9 runs at <= 2.9% (*Histoplasma ohiense*, ONT, 2.7% measured, k-mer survival 0.5 or below). `error_est` of the runs: see `long_reads_real.tsv`/`results.jsonl` (`qc.error_est`).
The measured error is a lower bound of the true one, so a run at 2.7% measured may well be near the 3.3% the flag targets.

## The flag `high_error_reads` (derived from the sweep)

* **Signal.** The abundance histogram of a read set contains the k-mers seen at least twice; an error k-mer is almost always seen once. So the peak (`coverage_est`, the base coverage of the error-free k-mers) is the true base coverage times the k-mer survival (1-e)^k. With the base coverage = total bases / size of the nearest reference
  genome (`n_hashes` x scaled; the size of a related genome is within about 20%, which moves the estimate by under 1 point of error): survival = coverage_est / base coverage, e = 1 - survival^(1/21). Sweep medians of the survival (`coverage_est` over the true depth, table at the end): 0.9 (0.5% error), 0.8 (1%), 0.6-0.75 (2%), 0.2-0.4 (5%); from 8% error the peak is
  hard to find at all (0.0-0.2), which also reads "high error".
* **Rule.** Reads only: survival <= 0.5 (error >= 3.3%) and base coverage >= 1x. Below 1x the estimate is noise (2% error at 0.25x: 58% of the inputs would be flagged, using the true genome size), so no flag is set there.
* **Measured with the implementation** (`scripts/hpcc/39_long_reads_flag_check.sbatch`: 6 of the 24 genomes re-simulated, own genome hidden, `long_reads_flag_check.tsv`; the number is flagged / 6 and the median `error_est`):

| error \ depth | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|
| 0.5% | 0/6 (0.7%) | 0/6 (0.8%) | 0/6 (1.4%) | 0/6 (0.9%) | 0/6 (0.8%) |
| 1% | 0/6 (0.8%) | 0/6 (1.4%) | 0/6 (1.8%) | 0/6 (1.4%) | 0/6 (1.4%) |
| 2% | 2/6 (2.0%) | 1/6 (2.1%) | 1/6 (2.7%) | 0/6 (2.7%) | 0/6 (2.4%) |
| 5% | 5/6 (7.5%) | 6/6 (5.0%) | 6/6 (5.7%) | 6/6 (6.9%) | 6/6 (5.2%) |
| 8% | 6/6 (14.6%) | 6/6 (15.8%) | 6/6 (8.2%) | 6/6 (8.8%) | 6/6 (16.3%) |
| 10% | 4/6 (18.3%) | 5/6 (18.7%) | 6/6 (10.2%) | 6/6 (10.7%) | 6/6 (14.6%) |
| 12% | 2/6 (20.0%) | 4/6 (21.2%) | 6/6 (11.3%) | 6/6 (12.2%) | 6/6 (14.5%) |
| 15% | 1/6 (24.4%) | 1/6 (24.7%) | 6/6 (25.6%) | 6/6 (14.4%) | 6/6 (16.4%) |


  **False positives: 0 of 60 read sets at 0.5% and 1% error (1x-30x); 4 of 30 at 2% (est. 2.0-2.7%, flagged where the true error is at the edge of the rule).** Short reads: none in the #32 simulations (survival 0.57-1.0 at 0.1-1% error; the flag also needs >= 1x), and no flag on 8 of the 9 real runs <= 2.9%.
  **Sensitivity:** from 5x every error >= 5% is flagged in 6/6 genomes; at 2x: 6, 6, 5, 4, 1 of 6 for 5, 8, 10, 12, 15% error; at 1x 5, 6, 4, 2, 1 of 6. Where the peak is missing (10-15% error with <= 2x) there is no `coverage_est`, so nothing is flagged: **no flag does not prove low error**, only a flag is informative. `error_est` is
  close up to 5% (2.0-2.7% for 2%, 5.0-7.5% for 5%) and overestimates above (8-16% for 8%, 10-19% for 10%, 14-26% for 15%): the histogram has no clear peak there, so the value says "very high error", not how high.
* **Unchanged:** every call, status and threshold; `results.tsv` has no new column (the flag is in `qc_flags`; `qc.error_est` and `qc.kmer_survival_est` are in `results.jsonl`).

## Why a smaller k or a finer scaled is not needed

At 15% error and 5x the genus is already called for 0.71 of the genomes (the plateau is 0.75-0.79) with precision 1.00, although a 21-mer survives with probability 0.03 (0.85^21) and each genome k-mer is then seen 0.17 times on average at 5x: the sketch keeps 1 hash in 10,000.
A smaller k (survival 0.85^15 = 0.09) was **not tested** (it needs a second index); it would not be needed for the genus, which is already at the plateau from 5x, and species calls at >= 8% error are limited by the species cutoff (ANI 0.99) more than by sensitivity: precision is 1.00 there and the calls that are missing are the ones a lower cutoff would add with lower precision.
The cost is a second index, so the decision is: do not build it. If long-read species calls at >= 8% error matter, correct the reads first (consensus or Dorado/`ccs` correction) and classify the corrected reads.

## Recommended use (copied to docs/limitations.md and docs/tutorial.md)

| reads | settings | what to expect |
|---|---|---|
| HiFi, ONT R10.4 (<= 3%) | defaults; `--min-abund auto` when the depth is >= 10x | genus from 0.5x, species from 1x with the usual 0.99 cutoff; precision 1.00 |
| ONT R9.4, 5-10% | defaults (min-abund 2); never `--min-abund 1` | genus from about 2-5x (10% error: 5x), species rarely; `high_error_reads` is set from >= 1x |
| PacBio CLR, 10-15% | defaults; at least 5x | genus only; the species and the call at < 2x are withheld |
| any | never `--min-abund 1` for long reads | ANI falls to 0.85-0.90 and every call is lost |

## Caveats

* The error model of the sweep is uniform and independent; real ONT/CLR errors cluster, so the real survival is a little higher than (1-e)^21 and the simulation is a conservative bound. The real runs (measured error 0.1-13%) are 20, from 20 genera; the platforms in the 4 measured classes are not balanced (CLR-like includes ONT 2019 and a Sequel II).
* The sweep and the runs use the full index; a user with only the reps index (`bfd-reps`, one genome per species) cannot get species calls for a species whose reps genome is hidden, so those validation numbers do not transfer to the reps index.
* Only fungi with chromosome-level assemblies; no metagenomic long reads (a mixture is not detected, #31).

## Tables of the sweep (24 genomes; own genome hidden unless stated)

**Genus called (share of the 24 genomes), own genome hidden, min-abund 2**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.62 | 0.75 | 0.71 | 0.75 | 0.79 | 0.79 | 0.75 |
| 1% | 0.71 | 0.75 | 0.79 | 0.75 | 0.75 | 0.75 | 0.75 |
| 2% | 0.58 | 0.79 | 0.79 | 0.79 | 0.79 | 0.75 | 0.75 |
| 5% | 0.46 | 0.58 | 0.75 | 0.75 | 0.75 | 0.75 | 0.71 |
| 8% | 0.21 | 0.50 | 0.71 | 0.71 | 0.71 | 0.75 | 0.71 |
| 10% | 0.12 | 0.29 | 0.54 | 0.67 | 0.79 | 0.75 | 0.71 |
| 12% | 0.04 | 0.21 | 0.38 | 0.58 | 0.79 | 0.75 | 0.75 |
| 15% | 0.00 | 0.04 | 0.21 | 0.42 | 0.71 | 0.75 | 0.75 |

**Genus precision of the calls made, own genome hidden, min-abund 2**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 5% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 8% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 10% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 12% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 15% | - | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |

**Phylum called, own genome hidden, min-abund 2**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.62 | 0.88 | 0.92 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1% | 0.75 | 0.92 | 0.92 | 0.96 | 1.00 | 1.00 | 1.00 |
| 2% | 0.58 | 0.88 | 0.92 | 0.96 | 1.00 | 1.00 | 1.00 |
| 5% | 0.46 | 0.58 | 0.79 | 0.92 | 1.00 | 1.00 | 1.00 |
| 8% | 0.21 | 0.50 | 0.71 | 0.79 | 0.92 | 1.00 | 1.00 |
| 10% | 0.12 | 0.29 | 0.54 | 0.75 | 0.92 | 0.96 | 0.96 |
| 12% | 0.04 | 0.21 | 0.38 | 0.58 | 0.92 | 0.92 | 1.00 |
| 15% | 0.00 | 0.04 | 0.21 | 0.42 | 0.79 | 0.88 | 1.00 |

**Species recall over the 17 genomes with a conspecific in the reference, min-abund 2 (default)**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 | 0.65 |
| 1% | 0.76 | 0.76 | 0.82 | 0.82 | 0.82 | 0.76 | 0.29 |
| 2% | 0.76 | 0.76 | 0.82 | 0.82 | 0.76 | 0.65 | 0.00 |
| 5% | 0.53 | 0.65 | 0.76 | 0.82 | 0.71 | 0.35 | 0.00 |
| 8% | 0.18 | 0.41 | 0.59 | 0.47 | 0.53 | 0.35 | 0.00 |
| 10% | 0.12 | 0.24 | 0.29 | 0.47 | 0.35 | 0.24 | 0.00 |
| 12% | 0.06 | 0.00 | 0.18 | 0.24 | 0.24 | 0.24 | 0.00 |
| 15% | 0.00 | 0.00 | 0.00 | 0.24 | 0.18 | 0.00 | 0.00 |

**Species recall, min-abund 1**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.65 | 0.65 | 0.59 | 0.06 | 0.00 | 0.00 | 0.00 |
| 1% | 0.06 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 2% | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 5% | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 8% | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 10% | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 12% | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 15% | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

**Species recall, min-abund auto**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 |
| 1% | 0.76 | 0.76 | 0.82 | 0.82 | 0.82 | 0.82 | 0.82 |
| 2% | 0.76 | 0.76 | 0.82 | 0.82 | 0.76 | 0.76 | 0.82 |
| 5% | 0.53 | 0.65 | 0.76 | 0.82 | 0.71 | 0.35 | 0.76 |
| 8% | 0.18 | 0.41 | 0.59 | 0.47 | 0.53 | 0.35 | 0.00 |
| 10% | 0.12 | 0.24 | 0.29 | 0.47 | 0.35 | 0.24 | 0.00 |
| 12% | 0.06 | 0.00 | 0.18 | 0.24 | 0.24 | 0.24 | 0.00 |
| 15% | 0.00 | 0.00 | 0.00 | 0.24 | 0.18 | 0.00 | 0.00 |

**Species precision of the calls made, min-abund 2**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 1% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | - |
| 5% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | - |
| 8% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | - |
| 10% | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | - |
| 12% | 1.00 | - | 1.00 | 1.00 | 1.00 | 1.00 | - |
| 15% | - | - | - | 1.00 | 1.00 | - | - |

**Genus called with the species hidden (genus precision follows), min-abund 2**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.38 | 0.46 | 0.42 | 0.46 | 0.50 | 0.50 | 0.46 |
| 1% | 0.42 | 0.46 | 0.50 | 0.46 | 0.46 | 0.46 | 0.46 |
| 2% | 0.33 | 0.50 | 0.50 | 0.50 | 0.50 | 0.46 | 0.46 |
| 5% | 0.33 | 0.33 | 0.46 | 0.46 | 0.46 | 0.46 | 0.38 |
| 8% | 0.08 | 0.25 | 0.38 | 0.42 | 0.42 | 0.46 | 0.38 |
| 10% | 0.00 | 0.17 | 0.33 | 0.42 | 0.50 | 0.46 | 0.42 |
| 12% | 0.00 | 0.08 | 0.25 | 0.33 | 0.54 | 0.46 | 0.46 |
| 15% | 0.00 | 0.00 | 0.08 | 0.25 | 0.38 | 0.46 | 0.42 |

**Genus precision with the species hidden**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.89 | 0.91 | 0.90 | 0.91 | 0.92 | 0.92 | 0.91 |
| 1% | 0.90 | 0.91 | 0.92 | 0.91 | 0.91 | 0.91 | 0.91 |
| 2% | 1.00 | 0.92 | 0.92 | 0.92 | 0.92 | 0.91 | 0.91 |
| 5% | 1.00 | 1.00 | 0.91 | 0.91 | 0.91 | 0.91 | 0.89 |
| 8% | 1.00 | 1.00 | 1.00 | 0.90 | 0.90 | 0.91 | 0.89 |
| 10% | - | 1.00 | 1.00 | 0.90 | 0.92 | 0.91 | 0.90 |
| 12% | - | 1.00 | 1.00 | 0.88 | 0.85 | 0.91 | 0.91 |
| 15% | - | - | 1.00 | 1.00 | 1.00 | 0.82 | 0.90 |

**Genus called with the genus hidden (every call is wrong), min-abund 2**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 |
| 1% | 0.08 | 0.04 | 0.08 | 0.08 | 0.08 | 0.08 | 0.04 |
| 2% | 0.04 | 0.08 | 0.08 | 0.08 | 0.08 | 0.08 | 0.04 |
| 5% | 0.00 | 0.00 | 0.04 | 0.04 | 0.08 | 0.04 | 0.04 |
| 8% | 0.00 | 0.00 | 0.00 | 0.04 | 0.04 | 0.04 | 0.04 |
| 10% | 0.00 | 0.00 | 0.00 | 0.04 | 0.04 | 0.04 | 0.04 |
| 12% | 0.00 | 0.00 | 0.00 | 0.04 | 0.08 | 0.04 | 0.04 |
| 15% | 0.00 | 0.00 | 0.00 | 0.00 | 0.08 | 0.12 | 0.04 |

**Median ANI of the top hit, own genome hidden, min-abund 2**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.996 | 0.995 | 0.995 | 0.996 | 0.996 | 0.995 | 0.991 |
| 1% | 0.995 | 0.994 | 0.995 | 0.996 | 0.996 | 0.994 | 0.985 |
| 2% | 0.996 | 0.994 | 0.995 | 0.995 | 0.995 | 0.991 | 0.974 |
| 5% | 0.997 | 0.995 | 0.994 | 0.992 | 0.991 | 0.987 | 0.961 |
| 8% | 0.989 | 0.990 | 0.989 | 0.989 | 0.986 | 0.984 | 0.963 |
| 10% | 1.000 | 0.992 | 0.988 | 0.985 | 0.983 | 0.981 | 0.967 |
| 12% | 1.000 | 0.988 | 0.984 | 0.979 | 0.980 | 0.977 | 0.967 |
| 15% | 1.000 | 1.000 | 0.981 | 0.972 | 0.971 | 0.971 | 0.965 |

**Median ANI, min-abund 1**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.990 | 0.990 | 0.989 | 0.985 | 0.976 | 0.962 | 0.931 |
| 1% | 0.986 | 0.985 | 0.983 | 0.977 | 0.963 | 0.945 | 0.909 |
| 2% | 0.974 | 0.973 | 0.970 | 0.963 | 0.947 | 0.926 | 0.887 |
| 5% | 0.944 | 0.941 | 0.939 | 0.935 | 0.920 | 0.901 | 0.862 |
| 8% | 0.914 | 0.913 | 0.912 | 0.910 | 0.900 | 0.887 | 0.851 |
| 10% | 0.896 | 0.899 | 0.893 | 0.892 | 0.887 | 0.876 | 0.846 |
| 12% | 0.879 | 0.879 | 0.877 | 0.876 | 0.870 | 0.864 | 0.842 |
| 15% | 0.844 | 0.848 | 0.848 | 0.848 | 0.847 | 0.843 | 0.828 |

**Median `coverage_est` over the true depth, min-abund 2**

| error \ depth | 0.25x | 0.5x | 1x | 2x | 5x | 10x | 30x |
|---|---|---|---|---|---|---|---|
| 0.5% | 0.96 | 0.84 | 0.96 | 0.91 | 0.89 | 0.90 | 0.90 |
| 1% | 0.76 | 0.77 | 0.82 | 0.83 | 0.77 | 0.80 | 0.80 |
| 2% | 0.46 | 0.74 | 0.74 | 0.71 | 0.61 | 0.60 | 0.67 |
| 5% | 0.16 | 0.12 | 0.37 | 0.37 | 0.34 | 0.22 | 0.33 |
| 8% | 0.08 | 0.06 | 0.04 | 0.20 | 0.18 | 0.14 | 0.02 |
| 10% | 0.08 | 0.02 | 0.02 | 0.13 | 0.12 | 0.10 | 0.04 |
| 12% | - | 0.02 | 0.01 | 0.01 | 0.08 | 0.07 | 0.04 |
| 15% | - | 0.02 | 0.01 | 0.01 | 0.00 | 0.04 | 0.03 |

