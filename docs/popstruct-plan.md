# Plan: `popstruct` (FASTQ set in, PCA, clusters and a first-pass tree out)

Status: draft, 2026-11-04. Issue #15. Written after the Coccidioides and *B. dendrobatidis* explorations in
`analysis/popgen_cocci/` (scripts 00-15). Not implemented.

## Goal

Give the CLI a set of FASTQ (or FASTA) files, optionally a label file, and get back: QC with flagged samples, a
distance matrix, PCA/PCoA, clusters, a first-pass tree, and a report. Labels give colors and let the computed
clusters be compared with what the user already believes (species, lineage, source).

## What exists and can be reused

- `gmg-core`: FracMinHash sketching with abundances, `.gsk` files, `auto_min_abund` (`coverage.rs`), containment and Jaccard
  helpers. Only dependency is a hash crate, so it stays wasm-friendly.
- `gmg run`: merges R1/R2 and lanes per sample, `--sample-sheet`, `--resume`, results.jsonl, html/md/pdf report.
- `gmg-report`: hand-written SVG charts (ANI histogram, heatmap, bars). No scatter plot, no tree drawing.
- Not in the workspace: linear algebra, tree building, clustering.

## Target interface

```
gmg popstruct reads/ [--sample-sheet S.tsv] --labels labels.tsv -o out/
              [--scaled 200] [--max-reads 6000000] [--min-cov 4] [--decontam] [--index bfd.gmgi]
              [--cut 0.01 | --k 3] [--tree nj|upgma|none] [--root ID]
              [--compare-tree T.nwk] [--compare-pcs scores.tsv] [--resume]
```

Outputs: `qc.tsv`, `ani.tsv`, `mash.tsv`, `pca.tsv`, `clusters.tsv`, `tree.nwk`, `label_check.tsv`, `popstruct.json`,
`report.{html,md,pdf}`.

## Stages (each from something measured)

1. **Sketch.** One sketch per sample (R1+R2 merged). Default scaled 200 (1000 derivable). Abundance filter by the valley
   below the main coverage peak. A first version that took the first local minimum picked noise inside the main peak for
   high-coverage samples and dropped about 27 of 559; the valley-below-peak rule fixed it.
2. **QC flags, not silent drops.** Low coverage (< 4x; 14 Cocci, 1 Bd samples). Estimated-size outliers (contaminated
   Cocci samples reached 57 Mb against about 27; 26 Bd samples sat at 11.8 Mb against about 21 and distorted k-mer PC1).
   Unexpected organism via the optional `--index` screen (CA10 and CA12 were Lecanicillium, NM_7898 Histoplasma). Flagged
   samples are excluded with a stated reason; `--keep-flagged` overrides.
3. **Decontamination (optional).** For flagged samples keep a hash only if it is in a target-genus reference or in at
   least 3 unflagged samples. Brought 31 Cocci samples from a median 32.6 Mb to 27.4 Mb.
4. **Distances.** Containment ANI (C^(1/k), C = shared / min size) and mash distance. Sorted-hash merges with rayon.
   Cost grows as n^2; first target is hundreds of samples, thousands need an inverted index or scaled 1000.
5. **Ordination.** Hash-presence PCA and PCoA of 1-ANI. Only the hash PCA has been compared with SNP PCA (top-axis
   principal-angle cosines 0.98-0.999 within Cocci species and for cleaned Bd). PCoA has not, so the default is decided
   after that comparison.
6. **Clusters.** UPGMA cut at the largest gap in merge heights, user `--cut`, plus k-means on top PCs. DAPC deferred.
7. **Tree.** NJ on mash distance, Newick, `--root`, midpoint rooting only after outliers are removed (it failed with a
   Histoplasma outlier present).
8. **Labels.** TSV with `id`, `group`, optional extra columns. Colors PCA and tree; scores clusters against labels (ARI,
   contingency table, purity); flags per-sample conflicts from median ANI to each label plus nearest-neighbor label (this
   found the two Coccidioides swaps); missing ids are "unlabelled".
9. **External comparison.** `--compare-tree` and `--compare-pcs`: exact-clade test per cluster, principal-angle cosines and
   Procrustes, distance Spearman.
10. **Resolution advisor.** A report panel stating when sketches cannot resolve within-group structure, quoting
    nearest-neighbor ANI and recommending a SNP method. Bd: nearest-neighbor ANI about 0.999, lineages 2/3/4 and 5/7
    indistinguishable, Spearman against tree distance 0.05 among the closest pairs.

## Report additions

Labeled scatter with hover ids; tree with colored tips; size-versus-coverage QC plot; heatmap downsampled or tree-ordered
beyond 60 samples; tables of flagged samples and label conflicts. Same report engine and formats as `run`.

## Build order (small PRs)

1. Run the existing Rust `auto_min_abund` (`coverage.rs`, first rise after abundance 2, then a peak test) on the abundance
   histograms of the 27 high-coverage Cocci samples that the first Python rule mishandled. If it also fails, port the
   valley-below-main-peak rule (bounded by peak/2) with a regression test. See Review notes: the Rust and Python rules differ
   and neither has been run against the other.
2. Core: QC, sparse distances, UPGMA, NJ, Newick; golden tests against the Python outputs on the 15-genome fixture.
3. Ordination: PCA and PCoA with a small eigensolver.
4. Labels: parser, ARI, conflict checks.
5. CLI subcommand, output schema doc, `--resume`.
6. Report charts.
7. External comparison and resolution advisor.
8. Validation: Cocci and Bd as documented HPCC regression checks; a simulated CI fixture (reads from a mutated genome
   under a coalescent tree, a few populations) because real data cannot ship; tutorial.

## Findings the design rests on

- Coccidioides: species gap in ANI (within about 0.995, between about 0.978), exact SNP clades for both species and for the
  19-sample immitis cluster; posadasii k=3 supported only for one split (12 of 262 samples discordant). Two `samples.csv`
  label problems found.
- Bd: nearly clonal, all pairs ANI >= 0.98; 5 coarse groups recovered (ARI 0.976, merged after seeing the ANI table, so not
  independent); cleaned k-mer PCA matches genome-wide VCF PCA (cosines 0.999, 0.995; distance Spearman 0.96).
- Distance Spearman against SNP-tree distance among close pairs is about 0.05-0.33: sketches give coarse structure only.

## Open decisions

1. Home: new crate `crates/popsketch` (lib + bin, depends on `gmg-core` and `gmg-report`, splittable later; recommended)
   versus a `gmg popstruct` subcommand now.
2. Paired reads: merge R1+R2 per sample (about 2x reading cost) or R1 only as in the tests.
3. PCA default: hash PCA (validated) or PCoA (cheaper, unvalidated).
4. Scale target: hundreds or thousands of samples.
5. Label file: one `group` column or several.

## Caveats

- All results are from two datasets, one clonal and one with clear species structure. The method is untested on
  recombining populations or on mixtures.
- Several choices were made after seeing the data (the 17 Mb size cut for Bd, the merged Bd groups, k for k-means).
- The SNP comparisons are not independent of the sketches: the SNPs came from mapped reads, likely the same read sets.

## Review notes (2026-11-04, Fable review; claims checked against the code)

Adopted:
- **PCA needs a different solver than PCoA.** PCoA is a dense symmetric eigendecomposition of the n x n matrix. Hash PCA can be
  done on the centred n x n Gram matrix built from the intersection counts already computed for distances (equivalent to the
  PCA when n is much smaller than the number of hashes), so no sparse SVD is needed. Linear algebra (nalgebra or faer, or a
  short Lanczos/subspace iteration) goes in the new crate behind an optional feature, never in `gmg-core`.
- **Crate layout.** Add the new crate to workspace `members` (the root `Cargo.toml` lists five crates and omits `gmg-py`).
  `rayon` lives in `gmg-cli`, not `gmg-core`. Keep UPGMA, NJ, Newick and ARI as plain code in a library with no heavy
  dependency. This supports the new-crate recommendation over a subcommand.
- **Scale.** Target a few hundred samples. About 150k hashes per 30 Mb genome at scaled 200 gives roughly 5e10 comparisons for
  n = 559, minutes with rayon. An inverted index is the right algorithm beyond that. Drop "thousands" as a target.
- **Resume covers sketching only.** `--resume` keys on file path and size (`run.rs`, around line 509). The abundance cutoff, QC
  flags, decontamination vote and everything after depend on the whole cohort and are recomputed every run.
- **Report engine reuse is smaller than stated.** The heatmap uses a fixed 34 px cell (`charts.rs`), unusable past about 60
  samples. Scatter plot, tree drawing and hover are new work.
- **Statistical guards to add.** (a) For each PC report its correlation with estimated size and coverage and warn when high
  (Bd k-mer PC1 was a size artifact). (b) Cluster cut: refuse to split unless the largest merge-height gap clearly exceeds the
  second largest, and report k = 1 otherwise, so a clonal population is not forced into two groups. (c) Mark decontaminated
  samples in every plot, exclude them from the resolution advisor, and keep decontamination optional: the vote removes private
  sequence by construction, so the "32.6 to 27.4 Mb" result measures size, not accuracy. (d) Treat the abundance cutoff as a QC
  output and expose it; a low-coverage sample falls back to 2.
- **Evidence still missing.** A negative control (shuffled labels, a single-population null), sensitivity to scaled, cutoff and
  sample size, and variance over repeated read subsamples. The cosines of 0.98-0.999 are for the top few axes and may mostly
  reflect the species or clade split; report the number of axes and variance share when quoting them. The two label swaps are
  supported by k-mer neighbours and, for 485B-1, by the SNP tree; B3476 is not in the SNP set, and neither is confirmed against
  the sample metadata.
- **Build order.** Move the simulated CI fixture to step 2 so every later step can be tested. Do the cutoff fix and the Gram-matrix
  PCA/PCoA together. Put the CLI subcommand before the report charts so something runs early. Defer the resolution advisor and
  external comparison.

Checked and qualified:
- The reviewer said the Rust cutoff has a known failure on noisy high-coverage histograms. I have not run it on those
  histograms; the Python rule that failed was my own first version (first local minimum, histogram truncated at abundance 59),
  and the Rust rule differs (it walks up from 2 until the smoothed count first rises, then requires a peak that clears
  `PEAK_MIN_HASHES` and `PEAK_OVER_VALLEY`). It may already handle those samples. Step 1 is now a test, not a port.
- The reviewer suggested reusing the read-set correction at `classify.rs:12`. That constant (`READS_SPECIES_MIN_ANI`) lowers a
  species cutoff for classification; it does not correct an ANI estimate, so it does not carry over to pairwise distances.
- The suggestion that min-size normalisation inflates the between-species ANI of 0.978 is a hypothesis. The two species have
  similar genome sizes, so I would test it before attributing the value to that.
