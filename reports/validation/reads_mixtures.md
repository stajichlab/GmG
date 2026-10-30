# Mixed read sets (#45, reads half): what GmG says about a two- or three-fungus sample today

Real fungal read pairs of two (or three) strains mixed in known shares, classified with the default read settings and with `--min-abund 1`. There is no mixture
detector yet (the decomposition is #31), so this measures what the single-lineage call, the top-hit ANI and the hit list reveal. Metrics that need a detector
(false-mixture rate on pure inputs, minor-component recall by a detector, abundance error, a detection limit of such a detector) are **not available**.
Index `bfd-2026.10.4-s10000.gmgi`, `known` thresholds, `gmg` from `main` at 7296a0da1 (built privately on the cluster; the shared binary predates `weak_phylum`).
Scripts: `scripts/mix_plan.py`, `scripts/mix_reads.py`, `scripts/mix_score.py`, `scripts/hpcc/32_mix_build.sbatch`, `33_mix_plan.sbatch`, `34_mixtures.sbatch`.
Per-input output stays on the HPCC (`GmG_runs/reads/mix`); all tables: [reads_mixtures_tables.md](reads_mixtures_tables.md).

## Design

* **Parents.** ZyGoLife paired Illumina read sets whose genus is called correctly alone at 1.5M pairs (reads_titration.md), plus *Coccidioides* read sets with a curated
  species. Truth = the lab's file-name genus resolved through the NCBI taxdump (as in reads_titration.md); species only for *Coccidioides*. 32 pairs, 8 per distance:
  `same_species` (two strains with the same species label; the labels of the ZyGoLife sets are not curated), `same_genus` (two named species), `same_order` (two genera),
  `diff_phylum`. 4 triples (60:30:10): two of three genera of one order, two of three phyla.
* **Shares.** Parent A at 100, 99.9, 99.5, 99, 95, 90, 80, 50 and 0% of the read pairs (0% = parent B alone; 100% = A alone: the pure controls). The 99.9, 99.5 and 99% rows
  are also the low-level cross-contamination test (0.1, 0.5, 1% of a second fungus in a pure sample).
* **Depth.** `low` = 600,000 read pairs (about 4x of a 45 Mb genome); `high` = 2,700,000 pairs (about 20x), for the 25 pairs whose parents both have that many.
  Each parent contributes the first pairs of its files. 517 inputs per setting (32 x 9 + 25 x 9 + 4); no errors in either run.
* Each cell has 5-8 pairs, so one pair is 12-20 points of a percentage. The ZyGoLife parents are early-diverging fungi (Mucoromycota, Chytridiomycota, Kickxellomycota),
  mostly thin clades in the reference: the results may differ for dense clades (*Aspergillus*, *Fusarium*).

## Results (default read settings: minimum k-mer abundance 2)

**The call at the first rank where the parents differ** (genus when they share it). Cells: the call is the larger parent's taxon / the smaller parent's taxon / no call, in % of
the pairs; at 50:50 either parent counts as correct. No input got a third taxon ("other" is 0 in every cell).

| A share | same_species, same_genus (genus) | same_order (family), low | same_order, high | diff_phylum (phylum), low | diff_phylum, high |
|---|---|---|---|---|---|
| 100% (A alone) | 100 / 0 / 0 | 100 / 0 / 0 | 100 / 0 / 0 | 100 / 0 / 0 | 100 / 0 / 0 |
| 99.9, 99.5, 99% | 100 / 0 / 0 | 100 / 0 / 0 | 100 / 0 / 0 | 100 / 0 / 0 | 100 / 0 / 0 |
| 95% | 100 / 0 / 0 | 100 / 0 / 0 | 100 / 0 / 0 | 88 / 12 / 0 | 80 / 20 / 0 |
| 90% | 100 / 0 / 0 | 100 / 0 / 0 | 71 / 0 / 29 | 88 / 12 / 0 | 80 / 20 / 0 |
| 80% | 100 / 0 / 0 | 75 / 12 / 12 | 43 / 29 / 29 | 75 / 25 / 0 | 80 / 20 / 0 |
| 50% | 100 (either) | 100 (either) | 100 (either) | 100 (either) | 100 (either) |
| 0% (B alone) | 100 | 100 | 100 | 100 | 100 |

* **Pure controls and cross-contamination.** A alone and B alone are called correctly at the decisive rank in every pair (and A's genus is called for every A alone), and a second fungus at 0.1, 0.5 or 1% does not change the call in
  any of the 32 pairs, at either depth. GmG is silent about contamination of this size (as designed: it has no detector); the only trace is the top-hit ANI (below).
* **Closely related parents do not disturb the call.** Two strains of a species, or two species of a genus, give the shared genus at every share: the minority only
  has the same k-mers as the majority. Nothing is lost, and nothing can be said about a second strain.
* **The minority wins in 12-25% of the cases at 95:5 to 80:20 when the parents are in different phyla**, and, for the same order, in 12% (4x) and 29% (20x) at 80:20. In the different-phylum
  table it is two pairs: *Apophysomyces* (Mucoromycota, thin in the reference) mixed with *Coccidioides* (Ascomycota, dense: many genomes) is called Ascomycota at
  95:5, 90:10 and 80:20 at both depths, although both parents alone have ANI 0.96-0.99; and *Cunninghamella* (pure ANI 0.92) with *Rhizoclosmatium* (0.97) is called Chytridiomycota
  at 80:20. The call follows the clade that has many agreeing reference genomes and the better ANI, not the larger read share. The minority never wins at 99:1 or beyond.
* **At 20x the same-order mixtures are called worse than at 4x** (71% instead of 100% at 90:10, 29% no call): with the default minimum abundance of 2 the minority
  is dropped when it is covered below about 1x (10% of 600,000 pairs of a 45 Mb genome is 0.4x, 10% of 2.7M is 2x), so at 4x it is invisible and does no harm, and at 20x it passes
  the filter and competes with the majority. The filter hides low-coverage foreign DNA, as in reads_contamination.md; the depth at which a minority becomes visible is
  set by the minority's coverage, not by its share.
* **No `not_fungal`, `distant_non_fungal_hit` or `weak_phylum` flag** on any of the 517 inputs (`weak_phylum` is not set for reads: sequencing errors lower their ANI).

**The top-hit ANI as a mixture meter.** For assemblies the ANI falls as p^(1/21) with the larger share p (stress_tests_2.md). For reads, ANI relative to the larger parent alone
(observed / predicted), default settings:

| A share | diff_phylum, low | diff_phylum, high | same_order, low | same_order, high | same_genus, low |
|---|---|---|---|---|---|
| 99% | 0.9971 / 0.9995 | 0.9958 / 0.9995 | 0.9990 / 0.9995 | 0.9984 / 0.9995 | 0.9994 / 0.9995 |
| 95% | 0.9955 / 0.9976 | 0.9933 / 0.9976 | 0.9993 / 0.9976 | 0.9929 / 0.9976 | 0.9981 / 0.9976 |
| 90% | 0.9937 / 0.9950 | 0.9912 / 0.9950 | 0.9993 / 0.9950 | 0.9916 / 0.9950 | 0.9987 / 0.9950 |
| 80% | 0.9889 / 0.9894 | 0.9862 / 0.9894 | 0.9996 / 0.9894 | 0.9955 / 0.9894 | 0.9998 / 0.9894 |

The law holds for parents of different phyla (observed slightly below predicted, by 0.0005-0.0043, at 80:20 to 99:1), but for reads the ANI loss is smaller than predicted whenever the parents are related:
the minority's reads still match the majority's reference genomes, so same-order and same-genus mixtures lose almost nothing. At 50:50 the best hit is one of the parents and the comparison with A alone is
meaningless. So an ANI below the species cutoff cannot tell a mixture from a divergent strain, and a clean ANI says nothing about a related minority.

**Is the minority taxon in the top-10 hits?** (share of pairs whose different-parent taxon at the decisive rank is among the first 10 reference hits; `hits.tsv`)

| minority share | same_order, low | same_order, high | diff_phylum, low | diff_phylum, high |
|---|---|---|---|---|
| 0.1% | 17% | 0% | 12% | 20% |
| 0.5% | 33% | 40% | 62% | 100% |
| 1% | 50% | 40% | 75% | 100% |
| 5% | 67% | 80% | 88% | 100% |
| 10% | 83% | 80% | 100% | 100% |
| 20% | 83% | 100% | 100% | 100% |
| 50% | 100% | 100% | 100% | 100% |

This is far better than in the assembly chimeras (reports/validation/stress_tests_2.md: 1-2 of 8 at 10-30%), but the reason is the parents: they are early-diverging fungi, thin clades
whose 10 best hits reach well beyond the genus, so a small minority can enter. It is not a property of reads. In the triples (below) the dense *Coccidioides* clade fills the top-10
list with *Coccidioides* genomes and the 30% parent is absent for both three-phylum triples. A hit list is not a mixture detector, and its sensitivity depends on how crowded the majority's clade is.

**Three genomes (60:30:10), default settings:** the call is the 60% parent's genus for all four triples; the 30% parent's genus is never the call. The 30% parent's genus is in the top-10 hits for the two
same-order triples and for neither three-phylum triple; the 10% parent's genus is in the top-10 hits for 1 of 4 (a same-order triple).

## Minimum abundance 1 versus 2

`--min-abund 1` (every k-mer counted) is worse on every measure and is not usable for reads (also reads_contamination.md):

* The pure controls lose calls: 12% of the same-species pure inputs at 4x and 14% of the same-genus pure inputs at 20x get no genus call (default: 0%).
* The minority wins more often: different-phylum at 90:10 and 80:20 and 4x: 50% of the pairs (default: 12 and 25%); at 20x: 20 and 40% (20, 20%); same-order at 20x, 80:20: 57% (29%).
* Three genomes: for 2 of the 4 triples the call is not the 60% parent's genus (one is the 30% parent's genus, one stops at the order).

## What the species call does here

The *Coccidioides* pairs (2 same-species, 1 same-genus) cannot test species retention: the species is called for only 5 of the 12 pure controls at the default setting (the 5 have ANI 0.991-0.997; the other 7 have ANI 0.957-0.985, below the
reads species cutoff of 0.99); the genus is called for all. Species behaviour in mixtures therefore rests on reads_contamination.md (assemblies: contamination.md): species is withheld first, never wrong.

## Decision

No change to the thresholds or the reads settings. Documented in docs/limitations.md: a read set from two related fungi is called to the shared taxon; a second fungus of a different phylum at 5-20% of the reads wins the phylum call when its clade is
dense in the reference and the other is thin (2 of 8 pairs here); one at 1% or less never changes the call and cannot be seen except as a small ANI loss; and the call alone cannot say that a sample is mixed. The decomposition (#31) is the remaining work for this issue.

## Not covered

Read sets with sequencing errors in the minority only (the minority here is real reads of the same quality), plants and animals as the second organism (reads_contamination.md covers bacteria and human), dense-clade parents on both sides (*Aspergillus* pairs; the
*Coccidioides* sets are the only dense clade here, and they are called to genus only), long reads, and RNA-seq mixtures.
