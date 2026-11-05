# popgen_cocci: exploratory scripts for issue #15

Python and shell prototypes used to test whether FracMinHash sketches of read sets (k=21, abundance tracked) can group
strains, give a PCA and a first-pass tree. Not part of the build and not tested in CI. Design: [docs/popstruct-plan.md](../../docs/popstruct-plan.md).

They were run on the UCR HPCC against two read sets (559 *Coccidioides* strains, 334 *B. dendrobatidis* strains) and
compared with existing SNP trees and VCFs. The paths in the scripts are those of that cluster and need editing to reuse
them. Raw data and results are not in this repository.

| step | script | what it does |
|---|---|---|
| 0 | `00_manifest.py` | sample sheet to a manifest of read pairs |
| 1 | `01_sketch.sbatch` | `gmg sketch --abund` at scaled 200 and 1000, one array task per sample |
| 2 | `02_popdist.py` | abundance cutoff, optional decontamination, containment ANI, groups, PCoA, hash PCA, NJ, label check |
| 2a-c | `02a_refhashes.py`, `02b_flags.py`, `02c_refsketch200.sh` | reference hashes for decontamination, flagged-sample list |
| 3 | `03_plots.py` | overview figure |
| 4-7 | `04_run_clean.sh`, `05_tree_check.py`, `06_guess_outliers.sh`, `07_final.sh` | clean runs, tree monophyly check, identify outliers with `gmg guess`, final per-species runs |
| 8 | `08_snp_concordance.py` | clusters and species against existing SNP trees |
| 9-10 | `09_bd_eval.py`, `10_bd_coarse.py` | Bd lineages: scoring and merged-group scoring |
| 11 | `11_snp_pca.py` | k-mer PCA against SNP PCA from alignments |
| 12-15 | `12_bdvcf_probe.sh`, `13a_bd_vcf_dump.sh`, `13b_bd_vcf_pca.py`, `14_bd_kmer_pc1.py`, `15_bd_rerun.sh` | Bd: genome-wide VCF PCA, what drives k-mer PC1, rerun without half-content samples |
| 16-17 | `16_pcoa_vs_snp.py`, `17_gap_rule.py` | PCoA vs hash PCA against SNP PCA; where a largest-gap UPGMA cut would fall |

**Known problem in the early scripts.** The built-in abundance cutoff in `02_popdist.py` returns 2 for almost every sample
(520 of 559 Coccidioides, all 334 Bd), which leaves error k-mers in and inflated sample sizes and distances. Use `--cutoffs`
with the output of a small probe of `gmg-core`'s `auto_min_abund` (see `18_cutoff_compare.py`; the probe is a scratch crate
that reads `.gsk` files and calls `gmg_core::coverage::auto_min_abund`). Results from before this was found are marked in
the plan; scripts 19-21 use the corrected cutoff.
| 18-21 | `18_cutoff_compare.py`, `19_rerun_rustcut.sh`, `20_splithalf_*.sh`, `21_splithalf_analyze.py` | Rust vs Python abundance cutoff; reruns with the Rust cutoff; split-half noise floor |

Findings, with their caveats, are summarised in the plan. Several choices (size cut, merged Bd groups, k for k-means)
were made after looking at the data.
