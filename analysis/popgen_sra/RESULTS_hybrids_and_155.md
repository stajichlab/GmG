# C. neoformans 1.4-1.8x group (#159) and P. oryzae host structure (#155), run 2026-11-17

## The *C. neoformans* runs at 1.4 to 1.8 times the assembly size (11 in the table; the 10 not already explained were mapped, SRR30534096 is *A. fumigatus*)
First 3M read pairs of each run mapped (minimap2 -x sr, `10_classify_changed.sbatch` with `CLASSDIR`) to H99 (*C. neoformans*, GCF_000149245.1) and JEC21 (*C. deneoformans*, GCF_000091045.1). Two clean runs of the same set (DRR829882, DRR829879) are the controls.

| Group | mapped to H99 | mapped to JEC21 |
|---|---|---|
| controls (2) | 97 to 98% | 58 to 60% |
| the 10 runs | 39 to 81% (8 of 10 between 64 and 81%) | 61 to 73% (9 of 10 between 65 and 73%; SRR34066256 18%) |

The group loses about 25 to 30 points on H99 compared with the controls and gains about 5 to 12 on JEC21: consistent with hybrids that carry an H99-like and a JEC21-like genome, which is also why they sit at about 1.4 to 1.8 times a single assembly. It does not prove it (a lineage divergent from both references would also lower both rates, and 3M pairs and a short-read mapper give coarse numbers); SRR34066256 (39% to H99, 18% to JEC21) is not explained by this and stays open. Results are on the HPCC under `GmG_popgen/classify_h99` and `classify_dene`. Hybrid genomes are a known case in *C. neoformans*, so `size_high` for them is a property of the sample, not an error of the estimator.

## #155: the P7 *P. oryzae* set against host (PopSeq catalog host_name, 281 runs, 158 with a host that is ok)
`popsketch run --labels host_labels.tsv --strata study` (labels from `15_p7_host_labels.py`; hosts with fewer than 5 runs pooled as `other_host`).
- Host explains much of the distances: PERMANOVA R2 0.68, p = 0.001 (999 permutations, shuffled within study); each host falls into a single cluster (rice 88/88, wheat and *Lolium* together 38/38, finger millet 8/8; ARI 0.80 at k = 7; the pooled `other_host` runs spread over seven clusters).
- The gap rule still finds no clear cut (jump 0.86 floors), so the structure is a set of lineages without a single sharp gap.
- Caveats: PERMDISP is also significant (groups differ in spread, so PERMANOVA alone is not evidence of different centres); only 37% of samples vary in host within a study, so the p-value rests on those; the 109 runs without a host and the `other_host` pool are not tested; a lineage and its host are confounded in this species.

What this does and does not settle for #155: the sketch recovers the known host-associated lineage structure of *P. oryzae* in a real set, with no fabricated clusters inside the rice group. It is not the test #155 asks for (a recombining population whose structure is known from SNPs): that still needs a SNP-based result for the same runs (`--compare-pcs` / `--compare-tree`).
