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
| `collapse.tsv` | with `--collapse`: every sample that entered the collapse, the representative that stands for it, the number of samples that representative stands for, the split-half distance used, and whether it is a representative |
| `split_half.tsv` | each sample's own split-half distance (1 - ANI between its two read halves), for the samples that carry halves: the technical noise of a comparison involving that sample (`--floor-samples` sets how many carry halves) |
| `group_tests.tsv` | for each label column: PERMANOVA (R2, pseudo-F, p), PERMDISP p, the stratum column, warnings |
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

**Do the labels explain the distances?** With a label file, `popsketch run` also tests each label column on the 1 - ANI matrix: PERMANOVA (does the label explain more of the squared distances than shuffled labels do, `--permutations`, default 999, 0 turns it off) and PERMDISP (do the groups differ in spread). `--strata COLUMN` names a label column, usually the study, and the permutations then shuffle labels only among samples with the same value there, because samples from one study are often near-clonal or sampled together and are not exchangeable across labels. Warnings: the label never varies inside a stratum (it is the study: no p-value); PERMANOVA and PERMDISP both significant (different spread, so not evidence of different centres); unequal group sizes with possibly unequal spread; groups under 5 samples; distances that are not Euclidean. Groups with one sample are left out. Results are in `group_tests.tsv`, `popsketch.json` (`group_tests`) and a report section; the page has fields for the strata column and the permutations. Limits: PERMDISP permutes the distances to the centroids by label (an approximation to permuting residuals); the test uses the samples in the analysis (flagged ones left out by default); it says nothing about causes. Code: `crates/popsketch/src/permanova.rs`; tests: its unit tests (calibration over 200 null sets, planted location and spread effects, within-study permutation) and `tests/fixture.rs`.

With a label file, the leading PCs of the hash PCA (`--dapc-pcs`, default 10; 0 turns it off) go into a discriminant analysis with the groups of the primary label column. It writes `dapc.tsv` (the discriminant scores, the group each sample is nearest to, and whether that is its own label), plots the first two axes, and lists the samples that are nearest to a group other than their label. This is a second hint, next to the ANI-based conflict check, that a label may be wrong; it needs at least 2 groups with 2 samples each and more samples than PCs.

## Collapsing near-identical samples: `--collapse`
Clonal or outbreak sets (for example most of the C. neoformans, P. oryzae and N. discreta sets) contain many samples that this sketch cannot tell apart, and they dominate the ordination, the tree and the label tests. `popsketch run --collapse [mean|max|median]` keeps one representative per group of near-identical samples and runs everything after the distances on the representatives (the clusters, PCA or PCoA, tree, DAPC and report). The label tests (`--labels`, PERMANOVA and PERMDISP) are the exception: they use all samples, because a collapsed cluster often mixes labels (in the clonal SRA sets 60-80% of the multi-sample clusters mix studies), so a representative's label would misstate its members. Two samples are near-identical when their distance (1 - ANI) is at most the technical noise of the comparison, which comes from each sample's own split-half distance (the distance between two sketches of alternate reads of the same sample): `mean` uses the mean of the two samples' distances (the default), `max` the larger, `median` the species median for every pair. A sample with no split halves gets the species median, so raise `--floor-samples` to give every sample halves (the default is 20). Samples are visited from the least noisy (best covered) to the noisiest and each unassigned one becomes a representative and takes the unassigned samples within the threshold of it; there is no chaining, and the result is one valid clustering (it depends on the visiting order, which is deterministic), not the smallest. `collapse.tsv` says which sample stands for which and the size of each cluster. (`permanova::test` can weight rows, for example by cluster size, and is tested for that; the pipeline does not use it because of the above.) `--collapse-within COLUMNS` (label columns, comma-separated, or `all`) makes the collapse label-aware: two samples are merged only when they carry the same labels in all of those columns (a missing label counts as its own value), so no representative stands for samples with different labels; the ordination and tree then show the labels honestly, at the cost of more representatives (see the plan, D1i). Without it, a cluster of near-identical samples from several studies is one point carrying one sample's label. If collapsing would leave fewer than 3 representatives it is skipped with a message. The thresholds were chosen from the split-half distances of seven SRA population sets (`bench/D1_collapse/`, `docs/plan_feasibility_tests.md` D1-D1e): collapse removes a lot only in clonal sets, and it is not a substitute for a cap on the number of samples. Not yet available in the web page. Code: `crates/popsketch/src/collapse.rs` (unit tests: representatives, noise order, no chaining, rules, determinism), `permanova.rs` (weights: weight 1 equals the unweighted test; a weight of k + 1 gives the same sums of squares as k exact copies), and `tests/fixture.rs`.

## Reading less: `--target-cov`

By default each sample is read up to `--max-bases`. `--target-cov C` instead stops reading a sample once the k-mer coverage estimated from the reads so far reaches C (checked every 200,000 reads, with R1 and R2 read in lockstep); `--max-bases` stays as a ceiling. A sample that never reaches the target is read up to the ceiling, as before, and gives the identical sketch. This helps when genome sizes in a set differ a lot, so that the same number of bases is very different coverage.

## Comparing with SNP results

`--compare-tree TREE.nwk` and `--compare-pcs SCORES.tsv` compare the run with a tree and PC scores you already have (for example from SNPs) for the samples in both (sample ids must match; characters Newick reserves are read as underscores). For the tree: the Spearman correlation of the tree's patristic distance with 1 - ANI (all pairs, and the closest half, where sketches are weakest), and for each cluster whether it is a clade of the tree and the best F1 against any clade or its complement. For the PCs: cosines of the principal angles between the top-2, top-3 and top-5 subspaces, the Procrustes disparity of the top 3, and the Spearman correlation of pairwise distances in the top-5 space. The PC table is tab- or comma-separated with an id column and numeric columns.

## In the browser

`crates/popsketch-wasm` builds the same analysis for WebAssembly (`pixi run popsketch-wasm`; about 1 MB), and `web/popsketch.html` is a page for it: choose FASTQ files or a folder and an optional label table, and get the report, the TSV and JSON files and the tree, with nothing uploaded. It differs from the command line in these ways:

- It reads only the first reads of each file (1,000,000 per file by default, set on the page), gunzipped with `DecompressionStream`, so the coverage and size flags depend on the genome size; read the coverage shown next to each flag.
- If your files are named by run accession (`SRR1713175_1.fastq.gz`), a button fetches country, study, host, strain and species for them from the ENA portal API and uses those as labels (only the accessions are sent). Reads are never fetched from an archive by the page: streaming the first reads of ENA FASTQs was tried, and ENA's file server was too slow and unreliable from a browser (a few KB/s with long waits in a handful of requests from one network, and AWS's copy sends no cross-origin header and holds `.sra` files popsketch cannot read), so archive fetching is a command-line job (`analysis/popgen_sra`).
- It runs on one thread (no `SharedArrayBuffer`, so no special server headers) and is meant for up to about 200 samples.
- No sketch cache or `--resume`, and no PDF. `--compare-tree` and `--compare-pcs` are there (the page takes a Newick file and a PC table). A `.gmgi` reference index can be loaded from a file on the page (the release index is about 213 MB, so the first load is slow) and adds the organism and mixture calls; without it the audit is cohort-relative only. The Markdown report is offered as one file with the figures embedded. The parity test checks the comparison, the Markdown text and figures, and the index screen against the command line.

The code is shared, not copied: `popsketch::pipeline::analyze` takes the sketches and the contents of the label file and returns the results with no file access; `run` is that plus discovery, sketching from files, the cache and the output files. The crate has a `parallel` feature (rayon) and a `native` feature (files, the command line, PDF), both on by default; the browser build turns both off. `popsketch_wasm::Session` sketches each sample from bytes and calls `analyze`. `crates/popsketch-wasm/tests/parity.mjs` (`pixi run popsketch-wasm-test`, also in CI) runs the simulated `popsketch demo` data through the native command and through the wasm build and requires the same sample counts, flags, noise floor, cluster verdict and partition, label check and ordination (numbers equal to 1e-9; the native code sums in parallel, the wasm code on one thread). The page was also driven in a headless browser with the demo reads: 32 samples, 27 in the analysis and a clear cut into 3 clusters, as on the command line.

## Sketching on its own, then analysing later

`popsketch sketch --id ID -o DIR R1.fastq.gz R2.fastq.gz` sketches one sample into a folder of sketches (including the split halves for the noise floor), and `popsketch run --from-sketches DIR -o out/` analyses every sample of the folder without reading any reads. The parameters (`-k`, `--scaled`, `--max-bases`) must match between the two; a mismatch is refused. This is for sets whose reads are not kept: one job per SRA run downloads, sketches and deletes the reads (`analysis/popgen_sra/03_fetch_sketch.sbatch`), and only the sketch stays: about 8 MB per sample at the default scale (the whole sketch about 3 MB, each split half about 2 MB; 2.7 GB for 300 *C. neoformans* samples, against hundreds of MB to several GB of reads each). The result is identical to analysing the reads directly (a test checks this).

### Sketch bundles for the web page

`scripts/popsketch_bundle.py SKETCH_DIR -o set.tar.gz [--labels labels.tsv] [--note TEXT]` packs a folder of sketches into one `.tar.gz`: a manifest (`popsketch-bundle.json`, format 1, with k, scale and the sample ids), the optional label table, and each sample's `<id>.json`, `<id>.gsk` and split-half sketches. It refuses a sample with a missing file, so a bundle is complete. The web page opens a bundle from your disk or from a link (`popsketch.html?data=URL`, the server must allow cross-origin requests) and does the audit, clusters, PCA, tree and report in the tab, so the reads are fetched and sketched once on a computer or cluster and the page is used to look at them and change thresholds. A bundle carries the split halves for the 20 samples the noise floor uses, not for all of them (the rest would not be used and are two thirds of the bytes), so it is about 3 MB per sample before gzip. The page reads sketches through `Session.add_sketch` (`sketch_from_bytes` in the crate), checked against `popsketch run --from-sketches` in the wasm parity test, which also found and fixed a difference: halves are used for the same pseudo-random subset of samples (`--floor-samples`) whoever made the sketches.

### How much of a run is enough

Is a sketch of part of a run the same as a sketch of the whole run? Tested on 12 runs each of *C. neoformans* (2.2 to 6.2 Gb of bases per run, whole-run k-mer coverage about 119x) and *Z. tritici* (about 37x), by sketching each run five ways (`analysis/popgen_sra/06_subsample_test.sbatch`, `07_subsample_compare.py`): the whole run, every Nth read pair to 1 Gb (`stride1g`, what the population sets use) and to 300 Mb, and the first 300 Mb and 100 Mb of reads. Distances (1 - ANI) of all 66 pairs against the whole run:

| | variant | Pearson r | mean abs. difference | nearest neighbour same | median coverage (whole run 119x / 37x) |
|---|---|---|---|---|---|
| *C. neoformans* | stride 1 Gb | 0.998 | 0.0081 | 75% | 31x |
| | stride 300 Mb | 0.998 | 0.0111 | 42% | 10x |
| | first 300 Mb | 0.998 | 0.0107 | 50% | 10x |
| | first 100 Mb | 0.998 | 0.0207 | 25% | 3x |
| *Z. tritici* | stride 1 Gb | 0.62 | 0.0016 | 33% | 13x |
| | stride 300 Mb | 0.09 | 0.0072 | 25% | 4x |
| | first 300 Mb | 0.21 | 0.0066 | 8% | 3x |
| | first 100 Mb | 0.20 | 0.0311 | 0% | 2x |

What this shows, and what it does not:

- Distances between lineages (*C. neoformans*, differences up to about 0.1) are preserved at every depth (r 0.998), but the numbers move: a 1 Gb sketch is off by about 0.008 on average, which is more than the noise floor (0.0025) of that species, so a distance from a subsample is not interchangeable with one from the whole run.
- Distances within a population (*Z. tritici*, mostly below 0.005) are mostly lost: at 1 Gb the correlation with the whole run is 0.62 and a sample's nearest neighbour is the same in one case in three; at 300 Mb or less it is essentially noise. A "no clear cut" for a recombining species from 1 Gb sketches therefore means the sketches cannot resolve fine structure, not that there is none. Use more of the run (not tested: 2 to 4 Gb) for within-population work. The depth series below qualifies this: in a second set of 12 deeper *Z. tritici* runs, whose distances span 0.01 to 0.04, sketches at 1, 2, 4 and 6 Gb agree with each other (r 0.998 or better), so the loss here is partly that this first set's distances were small compared with the sampling difference, and partly that the whole run is not a reliable reference when it is very deep.
- Taking the start of the reads (`head`) instead of every Nth pair made no clear difference at 300 Mb: both are limited by depth, not by position in the file (r 0.998 against 0.998 and 0.21 against 0.09 are within what 12 samples can show). This test cannot tell whether the first reads are biased at greater depth.
- The genome size estimate is stable down to 300 Mb (0.89 to 1.0 of the whole run's median) and breaks at 100 Mb (0.46 for *Z. tritici*).
- Limits: two species, 12 runs each, noise floors estimated from 12 samples, one subsample per run. It is a first look, not a validation.

**Depth ceiling (fixed in #157).** The coverage and genome size estimates come from the peak of the k-mer abundance histogram. Before the fix the peak test failed above about 100x k-mer coverage: at 150x and more the error k-mers seen 2 to 30 times (recurrent sequencing errors) outnumber the peak, which held only 37 to 46% of the hashes seen twice or more, below the 50% the test asked for. The sample then fell back to the low-coverage estimate (coverage about 0.6, size 40 to 200 Mb, abundance cutoff 2, so error k-mers stayed in the distances). In 12 *C. neoformans* runs (genome about 18 Mb) the peak was found for 11 at 1 Gb and 2 Gb, 5 at 4 Gb, 3 at 6 Gb and 1 for the whole run (about 400x). A peak is now also accepted when the hashes within 40% of its position hold at least 30% of the hashes seen twice or more (`PEAK_WINDOW_SHARE`); in the sets checked, false peaks hold 10% or less. With that, 10 of 12 *C. neoformans* runs have a peak at 4 Gb (was 5) and the mean distance difference between 4, 2 and 1 Gb falls from 0.011 to 0.001. On the 2,002 SRA sketches (1 Gb) four samples change from the low-coverage estimate to a peak and no cluster verdict or noise floor changes. Two of those get a size near 30 Mb (plausible), two get 3 to 4 Mb with a peak at 150x to 240x: most of their reads come from another small genome, and both were already flagged `unrelated_to_cohort`. Still wrong: a sample with a second genome (SRR27558955 estimates 45 to 120 Mb at every depth) and SRR10298347 at 4 Gb. The 6 Gb and whole-run sketches were deleted before the fix and have not been re-tested; the default `--max-bases 1000000000` is still the recommended depth, because it already gives 20 to 100x.

### Depth series: 1, 2, 4, 6 Gb and the whole run

Twelve deep runs each of *C. neoformans* (7 to 9 Gb, about 18 Mb genome) and *Z. tritici* (14 to 40 Gb, about 35 Mb genome) were each fetched once and sketched at 1, 2, 4 and 6 Gb of bases (every pair kept with the same probability, evenly spaced) and as the whole run (`08_depth_series.sbatch`, compared with `07_subsample_compare.py`; set `REF=s6g` to compare against another depth than the whole run). All five sketches of a run use the same parameters (`--max-bases 0`). Pairwise distances (1 - ANI) over the 66 pairs of 12 samples, with the 6 Gb sketch as the reference:

| | depth | Pearson r vs 6 Gb | mean abs. difference | nearest neighbour same | median k-mer coverage | runs with the abundance peak found |
|---|---|---|---|---|---|---|
| *Z. tritici* | 4 Gb | 0.9989 | 0.00011 | 92% | 74x | 12 of 12 |
| | 2 Gb | 0.9981 | 0.00015 | 83% | 37x | 12 of 12 |
| | 1 Gb | 0.9983 | 0.00019 | 83% | 18x | 12 of 12 |
| | whole run | -0.14 | 0.0196 | 8% | 0.6x | 2 of 12 |
| *C. neoformans* | 4 Gb | 0.78 | 0.0114 | 67% | 0.8x | 5 of 12 |
| | 2 Gb | 0.39 | 0.0236 | 33% | 71x | 11 of 12 |
| | 1 Gb | 0.42 | 0.0239 | 25% | 35x | 11 of 12 |
| | whole run | 0.88 | 0.0115 | 92% | 0.8x | 1 of 12 |

(The "peak found" count for the 6 Gb reference is 11 of 12 for *Z. tritici* and 3 of 12 for *C. neoformans*. The *Z. tritici* run SRR5829692 has no peak at 6 Gb and falls back to the low-coverage estimate.)

- ***Z. tritici***: from 1 Gb (about 18x) up, the distances do not change with depth beyond about 0.0002, which is below the noise floor of the sketches (0.0011 at 1 Gb), and the strain order is stable (83 to 92% same nearest neighbour). More bases did not buy resolution here. The whole run is the outlier, not the subsamples: at 15 to 40 Gb (400x or more) the peak is not found, the abundance cutoff falls to 2, error k-mers stay in, and every distance shifts by about 0.02. A whole run that deep is a worse answer than 1 Gb, not a better one.
- ***C. neoformans***: the table above was made before the estimator fix and does not show what depth does to resolution, because the estimate fails at 4 Gb and above and distances move with the cutoff. After the fix (see "Depth ceiling") 4, 2 and 1 Gb agree to a mean 0.001 (0.0011 and 0.0014 against 4 Gb), with two runs still estimated wrongly.
- So the practical guidance is the existing default: about 1 Gb (20 to 100x k-mer coverage for a fungal genome), not the whole run. The test did not show a gain from going deeper; it showed a loss for the very deep.
- Limits: two species, 12 runs each, one subsample per depth, and the 6 Gb reference is itself a subsample (and fails the estimate for 9 of 12 *C. neoformans* runs).

### Seven population sets from SRA

`analysis/popgen_sra/` fetched 2,002 runs (1 Gb of bases each, every Nth pair) for seven species and analysed each set (`04_run_species.sh`, `05_ordinations.sh`). The clusters are the gap rule at the default `--min-floors 3`; "country" and "study" are the labels from ENA (best adjusted Rand index over the candidate cuts reported by the label check). No clear cut means the largest gap between merge heights is not larger than a random tree of the same data (the null p-value is in `popsketch.json`).

| species | samples analysed | clusters | gap / noise floor (null p) | PC1, PC2 variance | best ARI: country, study |
|---|---|---|---|---|---|
| *Z. tritici* | 259 | no clear cut (1) | 0.4 (0.52) | 15%, 6% | 0.01, 0.14 |
| *C. neoformans* | 264 | clear: 2 (262 + 2) | 26.9 (0.06) | 64%, 18% | 0.01, 0.07 |
| *C. beticola* | 282 | no clear cut (1) | 0.1 (0.96) | 19%, 11% | 0.16, 0.07 |
| *P. nodorum* | 265 | no clear cut (1) | 0.4 (0.38) | 14%, 6% | 0.61, 0.22 |
| *B. cinerea* | 294 | no clear cut (1) | 0.6 (0.41) | 21%, 8% | 0.04, 0.10 |
| *N. discreta* | 203 | no clear cut (1) | 2.7 (0.78) | 43%, 27% | 0.19, 0.30 |
| *P. oryzae* | 251 | no clear cut (1) | 0.8 (0.76) | 64%, 13% | 0.27, 0.09 |

- Six of seven sets have no cut the gap rule trusts. For *Z. tritici* the depth series above says the 1 Gb sketches resolve distances of this size, so for that species "no clear cut" is a property of the population (a recombining one) and not of the sketch depth. It is not shown for the other species.
- *C. neoformans*: the one clear cut separates two samples, SRR34066268 and SRR34066269, both from PRJNA1279655, from the other 262. That is an outlier pair, not a lineage structure; the PC1 share (64%) is mostly those two. What they are has not been checked (a different species, a contaminant, or an aneuploid isolate are all possible). The label check reports 9 samples whose country label disagrees with their nearest neighbours (for example a Colombia-labelled sample closest to a China one) and 15 for study; these are listed in `label_check.tsv`.
- *P. nodorum*: country ARI 0.61 at 13 clusters may be mostly the sampling: 187 of 265 samples are Australian, so a country label here is close to a study label (study ARI is 0.22) and the two are not separated. Not checked beyond that.
- *P. oryzae*: PC1 holds 64% of the variance but the tree has no gap: a gradient, not groups. 5 country conflicts. 22 samples are flagged low coverage and 8 unrelated to the cohort (six are consecutive runs, ERR1976666 to ERR1976671).
- *N. discreta* has the largest gap relative to the noise floor (2.7) but a random tree of the same data gives a gap that large 78% of the time, so it is not a cut; its PCA shows two lobes (70% of the variance in two axes).
- Flags: low coverage is common (*P. nodorum* 34 of 265, *P. oryzae* 22, *Z. tritici* 39); `size_high` in *C. neoformans* (21) is at least in part the estimator failure described above for runs deeper than 100x and has not been checked sample by sample.
- Not checked: whether any of these sets contains duplicates, clones or mixed-species runs beyond what the flags say; the ordination figures and trees (`<slug>_pcoa/figures`, `tree.nwk`) were not inspected by eye for this note.

## Key options

`-k` (21), `-s/--scaled` (200), `--max-bases`, `--r1-only`, `--min-cov`, `--max-samples`, `--flagged`, `--index`, `--ordination`, `--components`, `--k`, `--cut`, `--min-floors`, `--no-tree`, `--root`, `--tree-floor`, `--floor-samples`, `--collapse`, `--collapse-within`, `--bitset-max-gib`, `--from-sketches`, `--target-cov`, `--size-z`, `--size-floor`, `--shoulder-z`, `--no-pca`, `--dapc-pcs`, `--strata`, `--permutations`, `--compare-tree`, `--compare-pcs`, `--report-format`, `--resume`, `--threads`. `popsketch demo` writes a simulated example. See `popsketch run --help`.

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

**A larger set.** Every fourth clean genome of the BFD set (2,500 genomes of many genera, no labels, no PCA; `analysis/popgen_cocci/37_popsketch_big_assemblies.sbatch`) ran to completion in 11 min 23 s with 24 threads and 19 GB of memory (about 3 times the time of the 1,598 *Aspergillus* genomes for 1.6 times the samples, because the genomes are larger and the all-pairs step grows with the square). 669 were excluded as unrelated to the cohort or for size, as expected for a set that is not one organism, and the rest were analysed. So the tool runs above the 1,500-sample design target and the 2,000 default limit (raised here with `--max-samples 3000`); the clusters of such a mixed set are not a meaningful result.

Not tested: a recombining population with known truth.
