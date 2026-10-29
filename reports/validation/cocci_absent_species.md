# Absent-species check for the reads species cutoff (#32)

Index built from BFD without the 6 *Coccidioides immitis* genomes; the 558 *Coccidioides* read sets classified against it (`scripts/hpcc/14_reads_absent_species.sbatch`). Reads of *C. immitis* (175 sets) must not be named *C. immitis*; the only relative is *C. posadasii*.

| reads | ANI of immitis reads to posadasii (q10 / median / q90) | immitis sets named at 0.9938 / 0.99 / 0.985 / 0.98 / 0.97 (all wrongly) |
|---|---|---|
| 500k | 0.9764 / 0.9783 / 0.9800 | 0.6% / 0.6% / 0.6% / 10.3% / 97.7% |
| 2M | 0.9764 / 0.9778 / 0.9788 | 0.6% / 0.6% / 0.6% / 0.6% / 97.7% |

True *C. posadasii* sets named at 0.99: 95.0% (500k), 95.6% (2M); at 0.9938: 91.9%, 85.4%. Limit: a single species pair; species complexes closer than ANI 0.98 would be named wrongly by this cutoff.
