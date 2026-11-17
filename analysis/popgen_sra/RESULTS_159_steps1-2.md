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

# Steps 3 and 4 (run 2026-11-16)

Scripts: `13_threshold_sweep.sbatch`, `14_threshold_eval.py` (step 3); `12_size_vs_assembly.py`, output `size_vs_assembly.tsv` (step 4). Ground truth for size is the species' assembly size in the PopSeq catalog; a call is right when est/assembly is 0.75 to 1.33.

## Step 3: the peak-window share (PEAK_WINDOW_SHARE, 0.30)
Every threshold from 0.20 to 0.40 was built from a copy of the source and run over all 2,002 finished sketches.

| Group | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 |
|---|---|---|---|---|---|
| all 2,002: right | 1942 | 1942 | 1941 | 1941 | 1941 |
| all: runs whose basis differs from 0.30 | 5 | 1 | 0 | 1 | 2 |
| tune (b_cinerea, c_beticola, c_neoformans): right /900 | 866 | 865 | 864 | 864 | 864 |
| held-out (n_discreta, p_nodorum, p_oryzae, z_tritici): right /1102 | 1076 | 1077 | 1077 | 1077 | 1077 |

The pick on the tuning species is 0.25, which scores the same on the held-out species as 0.30 (1077/1102). The estimate is insensitive to the threshold across 0.20 to 0.40: at most 5 of 2,002 runs change basis, and the share of right calls stays at 97.0%. **Keep 0.30**; the earlier worry that the choice rests on about 10 boundary samples is answered for these seven sets (the boundary samples are rare, and moving the threshold changes at most 5). Limits: the same seven fungal sets and one read type; ground truth is the assembly size, so a sample whose true size differs from the species' (hybrids, polyploids) counts as wrong here.

## Step 4: estimated size against assembly size, 2,002 runs with a size
62 runs (3.1%) are outside 0.75 to 1.33; **all 62 are already `excluded`**, so none enters an analysis. By species: b_cinerea 3, c_beticola 15, c_neoformans 18, n_discreta 3, p_nodorum 12, p_oryzae 8, z_tritici 2 (0.7% to 6.0%). The *C. neoformans* group at about 1.4 to 1.8 times the assembly size (26 to 34 Mb, 11 runs including ERR1520647 to ERR1520652) looks like a diploid or hybrid cluster, not contamination: not checked against the reads. Known cases in the list: SRR30534096 (mislabelled, 1.82; PopSeq stajichlab/PopSeq#52), SRR26462617 and ERR13824414 (non-fungal reads), SRR27558955 (second genome). Run `12_size_vs_assembly.py` for the full table.

## Still open in #159
Step 5 (SRR27558955 and SRR10298347 as separate cases; both grow with depth under the Poisson fallback, see the step 2 table), the unexplained 18 to 29 Mb step of DRR829884 between 1 and 2 Gb, and the *C. neoformans* 1.4 to 1.8x group.
