# #159 steps 1 and 2 (run 2026-11-16)

Scripts: `10_classify_changed.sbatch` (step 1), `11_deep_fixed.sbatch` (step 2). HPCC outputs: `/bigdata/stajichlab/jstajich/projects/GmG_popgen/{classify,depth,depth_cmp159}/` (purgeable). Fixed binary: `popsketch` 0.3.2.

## Step 1: reads of the four samples whose estimate changed (first 4M read pairs of each run)
| Run | Labelled | Mapped to the species' reference (minimap2 -x sr) | `gmg guess` (bfd_s1000) | Reading |
|---|---|---|---|---|
| ERR12084584 | *C. beticola* | 22% | *C. beticola*, 23.6k shared hashes | right species; about 78% of reads are another genome (likely host) |
| ERR13824414 | *Z. tritici* | 0.5% | no call (13 shared) | not *Z. tritici*; not in the fungal index |
| SRR26462617 | *C. neoformans* | 0.03% | no call (4 shared) | not *C. neoformans*; not in the fungal index |
| SRR30534096 | *C. neoformans* | 0.8% | *Aspergillus fumigatus*, 25.8k shared hashes | **mislabelled**; its ~30 Mb estimate fits *A. fumigatus* (~29 Mb), not *C. neoformans* (18.9 Mb) |

Supports the contaminant explanation for the 3-4 Mb estimates (ERR13824414, SRR26462617). Limits: the first 4M pairs only (a sorted run could bias this); the index holds fungi only, so what the non-fungal reads are is not identified.

## Step 2: deep runs re-sketched with the fixed binary (estimated size, Mb; basis)
| Run (reference size) | 1 Gb | 2 Gb | 4 Gb | 6 Gb | whole run |
|---|---|---|---|---|---|
| DRR829882 (*C. neoformans* 18.9) | 18.25 peak | 18.26 peak | 18.27 peak | 18.30 peak | 18.29 peak |
| SRR22188296 (*Z. tritici* 39.7) | 35.5 peak | 35.4 peak | 35.5 peak | 35.7 peak | 35.7 peak |
| SRR5829673 (*Z. tritici* 39.7) | 35.4 peak | 35.3 peak | 35.3 peak | 35.4 peak | 35.4 peak |
| DRR829884 (*C. neoformans*) | 18.3 peak | 28.9 peak | 29.3 peak | 29.4 peak | 29.3 peak |
| SRR10298347 (*C. neoformans*) | 17.9 peak | 17.9 peak | 48.5 poisson | 69.5 poisson | 96.4 poisson |
| SRR27558955 (*C. neoformans*) | 45.8 poisson | 73.1 poisson | 121.7 poisson | 164.7 poisson | 213.2 poisson |

- The fixed estimator is stable from 1 Gb to the whole run on the three clean samples, at 6 Gb and whole run (the cases deleted before the fix).
- DRR829884 jumps from 18 to 29 Mb between 1 and 2 Gb and stays: not explained (a second genome appearing at depth, or an estimator step); check its reads like step 1.
- The two known failures get worse with depth (error k-mers accumulate under the Poisson fallback): the estimate grows nearly linearly with bases. Step 5 of #159.
- Still open: steps 3 (held-out threshold), 4 (size vs assembly across 2,002 runs), 5.
