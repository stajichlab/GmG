# popsketch: audit a set of read files, then look at coarse structure

`popsketch run` takes a set of FASTQ files (or FASTA assemblies), optionally a label file, and writes:

1. **an audit**: which samples look wrong for the set (not the same organism, low coverage, a second genome mixed in, half the content missing, a label that the data contradict), each with the reason and the numbers;
2. **coarse structure** of the samples that pass: a hash-presence PCA, clusters, a first-pass tree and a report.

It compares k-mer sketches, not alignments. Sketches see large differences and miss small ones: use the tree to see which samples group together, not to say how strains are related. The report says where the data stop being informative (see "Resolution"). Design and evidence: [popstruct-plan.md](popstruct-plan.md).

```
popsketch run reads/ -o out/ --labels labels.tsv
popsketch run --sample-sheet samples.tsv -o out/ --flagged mark --index bfd.gmgi
```

## Inputs

- **Reads**: files or folders (searched recursively). The files of one sample (R1/R2, lanes) are merged by the same naming rules as `gmg run` (`--no-group` keeps them apart). A **sample sheet** (`--sample-sheet`, tab-separated, columns `id` and `path`; several rows with one id are one sample) names samples explicitly.
- **Assemblies** (FASTA) can be given too. They have no coverage, so the coverage and abundance checks do not apply to them; size and relatedness checks do.
- Each sample is read up to `--max-bases` (default 1 Gb, shared by its files; `0` = no limit). Equal bases give roughly equal coverage for similar genomes. `--r1-only` reads only the first file of each sample (fast).
- **Labels** (`--labels`): a table whose first column is the sample id and whose other columns are label sets (species, lineage, source ...). A blank or `NA` cell means unlabelled. The first label column is the primary one (it colours the plots); every column gets its own comparison.
- Up to `--max-samples` (default 1500) samples: all pairs are compared exactly, which grows with the square of the number of samples. Above the limit the run stops with a message instead of running for hours.

## What the audit flags

Every flag is relative to the cohort; nothing needs a reference. Each flag carries a reason such as `estimated size 11.8 Mb vs cohort median 21.2 Mb (z -9.6)`.

| flag | meaning |
|---|---|
| `low_coverage` | k-mer coverage below `--min-cov` (default 4), or too few repeated k-mers to estimate it. These samples lose k-mers by chance and look more distant than they are. |
| `size_low` / `size_high` | the estimated genome size (hashes above the abundance cutoff x scaled) is far from the cohort median (robust z above 5, with a floor of 3% of the median so a tight cohort does not flag everything). Low: content is missing. High: a second genome or extra content. |
| `abundance_shoulder` | k-mers at 0.1 to 0.6 times the main coverage exceed the cohort: a second organism at several-fold coverage. |
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
| `pca.tsv`, `pca_variance.txt` | PC scores (hash PCA by default; `--ordination pcoa` for PCoA of 1 - ANI), cluster per sample |
| `clusters.tsv` | cluster per sample (`NA` when there is no clear cut) |
| `tree.nwk` | neighbour-joining tree on mash distance, unrooted unless `--root ID`; branches shorter than the noise floor are collapsed |
| `label_check.tsv` | samples whose ANI matches another label better than their own |
| `popsketch.json` | everything above in one machine-readable summary |
| `sketches/` | the sketches; `--resume` reuses them (only sketching is resumed; QC, flags and everything after are recomputed) |

## Clusters

Samples are clustered by average linkage on 1 - ANI. The cut is made where the merge heights jump most, and is called **clear** only if that jump is at least `--min-floors` (default 3) times the noise floor. Otherwise the report lists the top candidate cuts and no cluster is assigned (use `--k N` or `--cut D`). If everything merges below the noise floor, the set is one population. The report also gives the probability of a gap ratio this large in a population with no structure: the largest-gap ratio alone is at least 2 in about 60% of structureless populations, so it is shown but is not the gate.

## Resolution

The **noise floor** is measured, not assumed: alternate reads of 20 samples go into two sketches and the distance (1 - ANI) between the halves is the technical noise; the 90th percentile is the floor (a fixed 0.001 is used if there are too few read sets to split). Differences smaller than the floor are not distinguishable from re-sampling the same reads: tree branches shorter than it are collapsed to polytomies, and the report states how many pairs are closer than the floor and whether each cluster's internal spread is above it. Where most pairs are inside the noise (near-clonal sets), the order of branches inside groups is not supported; use a SNP-based method for strain-level relationships.

Measured limits on two real sets: *Coccidioides* (species and a population cluster recover well, matching SNP trees on the top two axes) and *B. dendrobatidis* (nearly clonal: lineages 2, 3 and 4 are not distinguishable). Sketch distances among close pairs agreed poorly with SNP-tree distances (Spearman about 0.05 to 0.3).

## Key options

`-k` (21), `-s/--scaled` (200), `--max-bases`, `--r1-only`, `--min-cov`, `--max-samples`, `--flagged`, `--index`, `--ordination`, `--components`, `--k`, `--cut`, `--min-floors`, `--no-tree`, `--root`, `--tree-floor`, `--floor-samples`, `--resume`, `--threads`. See `popsketch run --help`.

## Validation on real read sets

Run on the 559 *Coccidioides* and 334 *B. dendrobatidis* read sets (R1 and R2 merged, 1 Gb per sample, release index) and compared with problems found earlier by other means (`analysis/popgen_cocci/29_check_popsketch.py`):

- ***Coccidioides***: the two species are one clear cut (k=2, a jump of 6.2 noise floors, adjusted Rand index 0.985 against the species label). The label check reports exactly the two samples already known to be mislabelled. All four non-*Coccidioides* read sets (*Histoplasma*, two *Lecanicillium*, *Fusarium*) are flagged, each by the cohort test as well as the index. All 20 mixtures that `gmg run` calls are flagged (20 of 20 with the index; the heavy ones by size or the shoulder alone). The 7% contaminant M223 and all lighter ones are found only through the index. The `center` label column, which does not follow the genetic structure, is reported as not checked instead of producing conflicts.
- ***B. dendrobatidis***: no clear cut, as expected for a near-clonal set (the top jump is 0.9 floors; the best candidate, 5 clusters, is the five coarse groups found by hand). 54% of all pairs are closer than the noise floor. Lineages 2, 3 and 4 are reported as indistinguishable (not checked). 26 of the 27 samples known to have half their content are flagged `size_low`, with 5 borderline extras; the index-confirmed 6.3% mixture is flagged. 12 samples are flagged `abundance_shoulder`; in a diploid this can be heterozygosity, not a second organism, and I have not separated the two.

Not tested: recombining populations, assemblies as input on real data, and sets of more than about 560 samples.
