# popsketch: audit a set of read files, then look at coarse structure

`popsketch run` takes a set of FASTQ files (or FASTA assemblies), optionally a label file, and writes:

1. **an audit**: which samples look wrong for the set (not the same organism, low coverage, a second genome mixed in, half the content missing, a label that the data contradict), each with the reason and the numbers;
2. **coarse structure** of the samples that pass: a hash-presence PCA, clusters, a first-pass tree and a report.

It compares k-mer sketches, not alignments. Sketches see large differences and miss small ones: use the tree to see which samples group together, not to say how strains are related. The report says where the data stop being informative (see "Resolution"). A worked example on a simulated set (`popsketch demo`) is in [popsketch-tutorial.md](popsketch-tutorial.md). Design and evidence: [popstruct-plan.md](popstruct-plan.md).

```
popsketch run reads/ -o out/ --labels labels.tsv
popsketch run --sample-sheet samples.tsv -o out/ --flagged mark --index bfd.gmgi
```

## Inputs

- **Reads**: files or folders (searched recursively). The files of one sample (R1/R2, lanes) are merged by the same naming rules as `gmg run` (`--no-group` keeps them apart). A **sample sheet** (`--sample-sheet`, tab-separated, columns `id` and `path`; several rows with one id are one sample) names samples explicitly.
- **Assemblies** (FASTA) can be given too. They have no coverage, so the coverage and abundance checks do not apply to them; size and relatedness checks do.
- Each sample is read up to `--max-bases` (default 1 Gb, shared by its files; `0` = no limit). Equal bases give roughly equal coverage for similar genomes. `--r1-only` reads only the first file of each sample (fast).
- **Labels** (`--labels`): a table whose first column is the sample id and whose other columns are label sets (species, lineage, source ...). A blank or `NA` cell means unlabelled. The first label column is the primary one (it colours the plots); every column gets its own comparison.
- Up to `--max-samples` (default 1500) samples: all pairs are compared exactly, which grows with the square of the number of samples. Above the limit the run stops with a message instead of running for hours. Measured on synthetic hash sets (`cargo test -p popsketch --release --test scale -- --ignored --nocapture`), after sketching: 1,500 samples take about 1 minute (distances 0.9 s, clustering 0.02 s, tree 1.3 s, hash PCA 30 s, PCoA 38 s) and 3,000 samples about 4.5 minutes (distances 2.3 s, tree 13 s, hash PCA and PCoA about 125 s each, dominated by the eigen-solver on data with no structure, its slowest case). Larger sets work with a higher `--max-samples`, slower; sketching time depends on the reads (about 45 s per sample per core for 1 Gb).

## What the audit flags

Every flag is relative to the cohort; nothing needs a reference. Each flag carries a reason such as `estimated size 11.8 Mb vs cohort median 21.2 Mb (z -9.6)`.

| flag | meaning |
|---|---|
| `low_coverage` | k-mer coverage below `--min-cov` (default 4), or too few repeated k-mers to estimate it. These samples lose k-mers by chance and look more distant than they are. |
| `size_low` / `size_high` | the estimated genome size (hashes above the abundance cutoff x scaled) is far from the cohort median (robust z above 5, with a floor of 3% of the median so a tight cohort does not flag everything). Low: content is missing. High: a second genome or extra content. |
| `abundance_shoulder` | k-mers at 0.1 to 0.4 times the main coverage exceed the cohort after removing the part that tracks the mass at 0.4 to 0.65 times (where a heterozygous peak sits): a second organism at several-fold coverage. A heterozygous or polyploid genome raises both and is mostly absorbed; a contaminant raises only the shoulder. If more than 2% of the cohort carries a shoulder the pattern points to heterozygosity or varying ploidy (as in *B. dendrobatidis*), so it is listed as a **note** and the samples are not excluded. A contaminant near half the coverage cannot be told from heterozygosity. Assessable from a k-mer coverage of about 12. |
| `unrelated_to_cohort` | the nearest other sample is far below the cohort's usual nearest-neighbour ANI: the sample shares almost nothing with the rest (another species or genus). |
| `too_few_hashes` | fewer than 500 hashes after the abundance cutoff. |
| `different_organism`, `mixture_suspected` | with `--index`: the index calls a different genus from most of the cohort, or finds a second taxon explaining at least 1% of the hashes. |

**What it cannot see.** Light contamination (roughly 1 to 7% of the reads) is invisible without a reference: a contaminant at 1x next to a 40x genome looks like sequencing error. On 558 *Coccidioides* read sets, 13 of 20 index-confirmed mixtures were at 1 to 3% and none of them was detected from the abundance spectrum or size; the heavy ones were (size for the 31 to 46% mixtures at low coverage, the shoulder for a 20% contaminant at 41x). `--index` adds the index's own mixture call for the light cases.

## What happens to a flagged sample (`--flagged`)

- `exclude` (default): left out of the distances, PCA, clusters and tree; still listed in `qc.tsv` and the report with its reason.
- `mark`: kept in everything and drawn with its own marker.
- `decontaminate`: for samples flagged only for contamination (`size_high`, `abundance_shoulder`, `mixture_suspected`), keep only hashes that at least 3 clean samples also carry, then put the sample back (drawn as a diamond). Other flagged samples stay out. This removes private sequence by construction and pulls a sample toward the cohort, so use it to look, not to measure.

## Outputs (`-o DIR`)

| file | content |
|---|---|
| `report.html` | one page: audit, plots, labels, structure, resolution statements |
| `qc.tsv` | every sample: status, flags, reasons, reads, bases, abundance cutoff, estimated size, coverage, index call |
| `ani.tsv`, `mash.tsv` | pairwise containment ANI and mash distance of the samples in the analysis |
| `dapc.tsv` | with a label file: DAPC scores on the primary label column, the group each sample is nearest to |
| `pca.tsv`, `pca_variance.txt` | PC scores (hash PCA by default; `--ordination pcoa` for PCoA of 1 - ANI), cluster per sample |
| `clusters.tsv` | cluster per sample (`NA` when there is no clear cut) |
| `tree.nwk` | neighbour-joining tree on mash distance, unrooted unless `--root ID`; branches shorter than the noise floor are collapsed |
| `label_check.tsv` | samples whose ANI matches another label better than their own |
| `popsketch.json` | everything above in one machine-readable summary |
| `report.md`, `figures/*.svg` | with `--report-format md` or `all`: the same report in Markdown with each figure as a standalone SVG |
| `report.pdf` | with `--report-format pdf` or `all`: the same report as a PDF (A4, selectable text, built with the pure-Rust renderer `gmg run` uses) |
| `sketches/` | the sketches; `--resume` reuses them (only sketching is resumed; QC, flags and everything after are recomputed) |

## Clusters

Samples are clustered by average linkage on 1 - ANI. The cut is made where the merge heights jump most, and is called **clear** only if that jump is at least `--min-floors` (default 3) times the noise floor. Otherwise the report lists the top candidate cuts and no cluster is assigned (use `--k N` or `--cut D`). If everything merges below the noise floor, the set is one population. The report also gives the probability of a gap ratio this large in a population with no structure: the largest-gap ratio alone is at least 2 in about 60% of structureless populations, so it is shown but is not the gate.

## Resolution

The **noise floor** is measured, not assumed: alternate reads of 20 samples go into two sketches and the distance (1 - ANI) between the halves is the technical noise; the 90th percentile is the floor (a fixed 0.001 is used if there are too few read sets to split). Differences smaller than the floor are not distinguishable from re-sampling the same reads: tree branches shorter than it are collapsed to polytomies, and the report states how many pairs are closer than the floor and whether each cluster's internal spread is above it. Where most pairs are inside the noise (near-clonal sets), the order of branches inside groups is not supported; use a SNP-based method for strain-level relationships.

Measured limits on two real sets: *Coccidioides* (species and a population cluster recover well, matching SNP trees on the top two axes) and *B. dendrobatidis* (nearly clonal: lineages 2, 3 and 4 are not distinguishable). Sketch distances among close pairs agreed poorly with SNP-tree distances (Spearman about 0.05 to 0.3).

## DAPC on the labels

With a label file, the leading PCs of the hash PCA (`--dapc-pcs`, default 10; 0 turns it off) go into a discriminant analysis with the groups of the primary label column. It writes `dapc.tsv` (the discriminant scores, the group each sample is nearest to, and whether that is its own label), plots the first two axes, and lists the samples that are nearest to a group other than their label. This is a second hint, next to the ANI-based conflict check, that a label may be wrong; it needs at least 2 groups with 2 samples each and more samples than PCs.

## Reading less: `--target-cov`

By default each sample is read up to `--max-bases`. `--target-cov C` instead stops reading a sample once the k-mer coverage estimated from the reads so far reaches C (checked every 200,000 reads, with R1 and R2 read in lockstep); `--max-bases` stays as a ceiling. A sample that never reaches the target is read up to the ceiling, as before, and gives the identical sketch. This helps when genome sizes in a set differ a lot, so that the same number of bases is very different coverage.

## Comparing with SNP results

`--compare-tree TREE.nwk` and `--compare-pcs SCORES.tsv` compare the run with a tree and PC scores you already have (for example from SNPs) for the samples in both (sample ids must match; characters Newick reserves are read as underscores). For the tree: the Spearman correlation of the tree's patristic distance with 1 - ANI (all pairs, and the closest half, where sketches are weakest), and for each cluster whether it is a clade of the tree and the best F1 against any clade or its complement. For the PCs: cosines of the principal angles between the top-2, top-3 and top-5 subspaces, the Procrustes disparity of the top 3, and the Spearman correlation of pairwise distances in the top-5 space. The PC table is tab- or comma-separated with an id column and numeric columns.

## Key options

`-k` (21), `-s/--scaled` (200), `--max-bases`, `--r1-only`, `--min-cov`, `--max-samples`, `--flagged`, `--index`, `--ordination`, `--components`, `--k`, `--cut`, `--min-floors`, `--no-tree`, `--root`, `--tree-floor`, `--floor-samples`, `--target-cov`, `--size-z`, `--size-floor`, `--shoulder-z`, `--no-pca`, `--dapc-pcs`, `--compare-tree`, `--compare-pcs`, `--report-format`, `--resume`, `--threads`. `popsketch demo` writes a simulated example. See `popsketch run --help`.

## Validation on real read sets

Run on the 559 *Coccidioides* and 334 *B. dendrobatidis* read sets (R1 and R2 merged, 1 Gb per sample, release index) and compared with problems found earlier by other means (`analysis/popgen_cocci/29_check_popsketch.py`):

- ***Coccidioides***: the two species are one clear cut (k=2, a jump of 6.2 noise floors, adjusted Rand index 0.985 against the species label). The label check reports exactly the two samples already known to be mislabelled. All four non-*Coccidioides* read sets (*Histoplasma*, two *Lecanicillium*, *Fusarium*) are flagged, each by the cohort test as well as the index. All 20 mixtures that `gmg run` calls are flagged (20 of 20 with the index; the heavy ones by size or the shoulder alone). M223 (a 7% contaminant) is flagged by the abundance shoulder as well as by the index; lighter ones (1 to 3%) are found only through the index. The `center` label column, which does not follow the genetic structure, is reported as not checked instead of producing conflicts.
- ***B. dendrobatidis***: no clear cut, as expected for a near-clonal set (the top jump is 0.9 floors; the best candidate, 5 clusters, is the five coarse groups found by hand). 54% of all pairs are closer than the noise floor. Lineages 2, 3 and 4 are reported as indistinguishable (not checked). 26 of the 27 samples known to have half their content are flagged `size_low`, with 5 borderline extras; the index-confirmed 6.3% mixture is flagged. the first version flagged 12 samples `abundance_shoulder` (the shoulder sat at 0.57 times the coverage in flagged and unflagged samples alike, and the flagged ones were barely more heterozygous in the VCF). After removing the part of the shoulder that tracks the heterozygous window, 5 remain. The index finds no mixture in them. Their allele balance at heterozygous SNP calls in the VCF (`analysis/popgen_cocci/32_bd_ploidy_from_ad.py`, against 12 unflagged controls) says that varying ploidy is not the explanation: three of the five (JEL412, Campana, TER_4.2) look like ordinary diploids (0.43 to 0.67 of their heterozygous calls near 0.5, a low-fraction tail of 6 to 7%, the same as the controls), so their shoulder is unexplained; two (BdBE3, KRBOOR_323) have about three times the heterozygous sites of the controls and a heavy low-fraction tail (48% and 30% against 7%), which fits a mixture of *B. dendrobatidis* genotypes or a lineage that maps poorly to the reference. One unflagged control (CCB1) looks triploid, so the shoulder does not track ploidy either way. Treat the shoulder in a diploid or polyploid organism as a hint to look at the sample, not as a finding.

The comparison options and DAPC were run on the same sets against SNP results made independently of the sketches (from mapped reads):

- ***Coccidioides*** against the 446-tip SNP tree (395 samples in both): the Spearman correlation of the tree's patristic distance with 1 - ANI is 0.85 over all pairs and 0.67 for the closest half of the pairs, and both clusters (the two species) are exact clades of the tree. Against the SNP PCA (393 samples in both): top-2 principal-angle cosines 0.999 and 0.989, Procrustes disparity of the top 3 of 0.016, distance correlation in the top-5 space 0.97. DAPC on the species label assigns 99.6% of the samples to their own species; the two it assigns to the other species are the two the label check reports.
- ***B. dendrobatidis*** against the genome-wide VCF PCA (241 samples in both): top-3 cosines 0.999, 0.996 and 0.978, Procrustes disparity 0.011, distance correlation 0.97 (the same numbers the earlier Python analysis gave). DAPC leaves out the lineage labels that the sketches cannot tell apart (lineages 2, 3 and 4) instead of reporting them as errors.

**A recombining population** (simulated, `tests/fixture.rs`): 24 samples that are mosaics of three founders in 4 kb blocks give no clear cut (the top jump is 0.7 floors), so recombination does not produce spurious clusters. Real within-species sets behave the same way (cached sketches, one species at a time): the 383 *C. posadasii* give no clear cut (top jump 1.9 floors, noise floor 0.0065), while the 175 *C. immitis* give one clear cut into 2 clusters (7.9 floors, noise floor 0.0019), as expected for a population with a deep split. These are partly recombining natural populations, not a controlled test, so the simulation remains the only case where the truth is known.

**Assemblies as input.** Run on genomes from the BFD clean-genome set (`analysis/popgen_cocci/33_make_assembly_inputs.py`, `34_popsketch_assemblies.sbatch`):

- 27 genomes (22 *Coccidioides*, 3 *Histoplasma* and 2 *A. fumigatus* as planted outgroups) take 3 s. Genus labels are recovered exactly (adjusted Rand index 1.0 at 3 clusters) and species labels exactly at 4. The planted outgroups are **not** flagged as samples: they form groups of their own, and the cohort-relative checks assume most of the set is one organism (the report warns that the set spans distant genomes). One *C. immitis* assembly (`GCA_982397225.1_27-FMI`, 870 contigs) is excluded as `unrelated_to_cohort` (nearest neighbour ANI 0.73 against a cohort median of 0.997); it has a normal size (24.9 Mb) and shares almost nothing with any other genome in the set. Checked with `gmg guess` against the release index: the only *Coccidioides* genome it resembles is itself (it is in the index under that label); its other hit is a sibling assembly of the same submission, `GCA_982397405.1_T55-F`, labelled *Ophidiomyces ophidiicola* (836 of 2,437 sketch hashes shared, which is roughly 95% ANI), and no other *C. immitis* genome appears among its hits. The submission is a mixed collection (a neighbouring accession, `4-FMI`, is a yeast, *Williopsis*). So the exclusion is a true positive: the *C. immitis* label of this assembly is very likely wrong, and it is the index's label that is in error here.
- 1,598 *Aspergillus* genomes (273 species) take 3 min 45 s with 24 threads and 14 GB of memory, so the 1,500-sample design target holds for assemblies. 3 empty files are reported as errors. As a mixed-species cohort it is outside the intended use: 84 divergent species are excluded as `unrelated_to_cohort` and two as `size_high`, the species label is recovered only coarsely (adjusted Rand index 0.66 at the 24 clusters of the largest gap, 0.76 at 195), and the label check proposes relabels that look like taxonomy corrections (for example *A. niger* to *A. costaricensis* and *A. flavus* to *A. oryzae*). The report now says so: when the deepest merge is at 1 - ANI of 0.10 or more (ANI below 90%) it adds a warning that the set spans distant genomes and that the cohort-relative checks assume one main organism. Use `popsketch` on one species or genus at a time; the exclusion of distant members is the cohort-relative check doing what it says.

**Threshold sensitivity** (`analysis/popgen_cocci/35_threshold_sweep.sbatch`, `36_summarise_sweep.py`; cached sketches, one threshold changed at a time):

| setting (default) | *Coccidioides*: flagged samples (34 by default) | *B. dendrobatidis*: flagged samples (43 by default) |
|---|---|---|
| `--min-floors` 1 to 5 (3) | 34, verdict `clear:2` throughout; 8 gives `no_clear_cut` | 43; `no_clear_cut` except 1, which accepts 5 clusters |
| `--size-z` 3 / 4 / 6 / 8 / 12 (5) | 39 / 34 / 33 / 33 / 30 | 44 / 43 / 36 / 33 / 33 |
| `--size-floor` 0.01 / 0.02 / 0.05 / 0.1 (0.03) | 48 / 37 / 33 / 30 | 51 / 43 / 33 / 30 |
| `--shoulder-z` 3 / 4 / 8 / 12 (5) | 30 / 36 / 32 / 31 | 38 / 43 / 42 / 40 |
| `--min-cov` 3 / 5 / 8 / 12 (5) | 31 / 34 / 40 / 56 | 43 / 43 / 43 / 44 |

The clustering verdict does not move with the QC thresholds, and the number of clusters stays at 2 for *Coccidioides* across the whole range. Flag counts move by a handful of samples for a factor of about 1.5 in `--size-z`, `--size-floor` or `--shoulder-z`, and more at the extremes: `--size-floor` 0.01 adds 14 to 8 size exclusions of borderline samples, and `--size-z` of 8 or more loses the 10 borderline *B. dendrobatidis* `size_low` calls. `--min-cov` 12 excludes 22 more *Coccidioides* because it drops everything sequenced below 12x. The `--min-floors` gate is what keeps the near-clonal *B. dendrobatidis* set from being cut: at 1 it reports 5 clusters although the top jump is only 1.5 floors. The defaults sit in a flat part of every curve for the species-level conclusions; the number of borderline flagged samples is the part that depends on them.

Not tested: a recombining population with known truth, and sets of more than about 1,600 samples.
