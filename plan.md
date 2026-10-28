# GmG (Guess My Genome): plan (draft for feasibility review)

Repo: GmG (renamed from PopSketch, 2026-10-27). **PopSketch** is now the working name for the later population-structure layer (strain trees, PCA) inside this repo; whether it stays here or splits out is decided after GmG works.

Status: DRAFT 2026-10-27. Numbers marked **[measured]** come from the 15-genome *A. fumigatus* test in `data/`. Everything marked **[est]** is an estimate that Milestone 0 must verify. Sections 12-14 record what was found on the HPCC and what is built (the `gmg-core` Rust crate); the rest is still plan.

## 1. Tools in scope

| Tool | Input | Output | Where it runs | Priority |
|---|---|---|---|---|
| **GmG: Guess My Genome** | assembled genome (FASTA) | species / species-group call, nearest references, containment/ANI, confidence | browser (wasm) + CLI | 1 |
| **GmG-reads** | FASTQ (unassembled) | same output, from reads | CLI/HPCC first; browser only for a read subsample | 2 |
| **PopSketch trees** | set of assemblies (later reads) | distance matrix, dendrogram, clusters | CLI/HPCC; viewer in browser | 3 |
| **MATPredict hook** | any of the above | group label that selects the MAT reference set | library/CLI call | consumer |

GmG and GmG-reads share one reference index and one scoring and threshold layer. PopSketch shares the sketching layer but a different comparison layer (all-vs-all vs. query-vs-DB). Keep those boundaries in the code from the start.

## 2. Goal 1 (GmG): refined aim

Original: "quickly guess what species group an unknown assembled genome falls into."

Refined, with testable criteria:

> Given a fungal assembly, return in under ~30 s in a browser (target; **[est]**) the best-matching species group from the BFD-scale reference set, with a calibrated confidence and an explicit "no confident match / novel" outcome.

Acceptance criteria (to be agreed before building):
- **Accuracy:** >= X% correct species call on held-out genomes (leave-one-out and leave-one-species-out). Proposed X = 95% at species level and 99% at genus level, to be revised after the first benchmark.
- **Novelty handling:** when the true species is absent from the DB, GmG must say so (report genus-level or "unknown" with the best ANI) and not force a species label.
- **Size budget:** the DB the browser must fetch is <= ~100 MB first load (cached afterwards). A larger optional tier is fetched lazily.
- **Latency:** sketch + search for one 30–40 Mb genome in wasm within the budget above. **[est]** Native sourmash sketching was ~1.2 s/genome here **[measured]**. wasm is typically 2–5x slower, so a few seconds is plausible but unverified.

Things to decide (these change the design):
1. What is a "species group"? NCBI species, a *sensu lato* complex (e.g., *A. fumigatus* complex), or a custom ANI-defined cluster (e.g., 95% ANI)? **Recommend:** ship both an NCBI-taxonomy label and an ANI-cluster ID. NCBI fungal taxonomy is inconsistent, and clusters are what the sketches can actually resolve.
2. Is the reference set the "BFD" ~20k genomes as-is, or dereplicated? Dereplication (e.g., one representative per ≥99% ANI cluster, with members kept in a second tier) is the main lever for DB size.
3. Taxonomy source for labels: NCBI taxdump via `sourmash tax` lineage CSVs. GTDB does not cover fungi.

## 3. Reference index design (the core engineering problem)

### What we know
- Sketch size scales as (genome length / scaled) × bytes per hash. At `k=21, scaled=1000`, the 15 test genomes (28–30 Mb) came to **4.1 MB total, ~270 KB each [measured]**.
- Estimates for 20k genomes at a ~35 Mb mean **[est]**:

| Setting | Per genome | 20k genomes |
|---|---|---|
| scaled=1000, 64-bit | ~280 KB | ~5.6 GB |
| scaled=1000, delta-coded | ~175 KB | ~3.5 GB |
| scaled=10000, 64-bit | ~28 KB | ~560 MB |
| scaled=10000, delta-coded | ~18 KB | ~350 MB |
| scaled=10000, 32-bit truncated | ~14 KB | ~280 MB |

- Species-level calls need much less resolution than strain-level. Strain-level resolution (99.2–99.99% ANI in the test set) needs scaled <= 1000.

### Proposed tiers

| Tier | Contents | Approx. size | Delivery | Used for |
|---|---|---|---|---|
| T0 | one representative per species/ANI-cluster, scaled=10000 | tens of MB **[est]** | static file, fetched on load | default GmG call |
| T1 | all genomes, scaled=10000, delta-coded | ~350 MB **[est]** | Releases/object store, HTTP range requests or opt-in download | within-genus refinement, CLI default |
| T2 | per-cluster scaled=1000 packs | GBs total, ~MBs per cluster | range-fetched per cluster | strain-level placement (PopSketch) |

### Format
- **CLI/HPCC:** native sourmash zip databases (compatible with `sourmash gather/search` and the branchwater plugin). Use an index/RocksDB-style inverted index for large searches. Check what the current plugin version supports.
- **Browser:** a small custom binary format (header, manifest JSON, then sorted delta-coded hash blocks per genome with an offset table) so that range requests can fetch single genomes or clusters. Build it from the sourmash zip with a script, so there is one source of truth.
- **Hash width:** 32-bit truncation is worth testing but adds false shared hashes (~0.3 per pair at scaled=1000 **[est]**). Gate it on benchmark accuracy.

### Hosting
- GitHub Pages is soft-limited to ~1 GB per site, and git rejects files over 100 MB. Put the app, the T0 index and the manifest on Pages. Put T1/T2 in Releases (2 GB/asset) or an object store (Zenodo, R2/S3, Hugging Face). **Verify CORS and range-request support on the chosen host before designing around it.**

### Build (HPCC)
- All genomes are already on the HPCC, so there is no download bottleneck.
- Sketching: ~1.2 s/genome **[measured]** gives ~7 core-hours for 20k genomes (one SLURM array job, scaled=1000 and 10000 in one pass via multiple `-p` params).
- Dereplication and clustering: `skani triangle` or branchwater `pairwise`/`cluster` at scaled=1000. Estimate minutes to hours. Do not hold the full 20k×20k matrix in the browser. Ship tree, cluster IDs and neighbor lists.
- Record provenance (accessions, assembly versions, tool versions, parameters) in a manifest shipped with each index release.

## 4. Search method (GmG scoring)

Do not only report one number. For each query:
1. Compute containment of query in each reference (and max containment). Containment is better than Jaccard here because assemblies differ in completeness and the reference may be larger.
2. Convert to containment-ANI (sourmash `--ani`) and rank.
3. Aggregate to taxon: best hit species, plus a vote over the top-N hits (guards against mislabeled references).
4. Report confidence from calibrated thresholds, not raw hit counts. Fit thresholds on the leave-one-out benchmark, per clade (thresholds differ between, say, *Aspergillus* and *Candida*).
5. Flag: low-ANI to everything (novel/unclassified), near-tie between species (complex), hits split across kingdoms (contamination: bacterial or human reads in the assembly).

Prior art to compare against before building anything novel (**check current versions and fungal coverage**): `sourmash gather/prefetch` and the branchwater plugin, `mash screen`, `skani search`, `fastANI`, `sylph` (low-coverage reads vs. genome DB with ANI), Kraken2/Bracken-style k-mer classifiers, and NCBI's own ANI-based taxonomy checks. The unique contribution of GmG is the fungal reference set plus the zero-install browser delivery, not a new algorithm. Say so in docs.

## 5. Goal 2 (GmG-reads): reads to genome guess

Reads change the problem:
- Sequencing errors make most unique k-mers singletons. **Abundance tracking and a minimum-abundance filter are required** (e.g., drop hashes seen once, or use a coverage-aware method).
- Coverage varies. At low coverage, containment of reads in a genome is small even for a perfect match. Use `gather`-style (greedy, containment-based) or `sylph`-style coverage-adjusted ANI, not the assembly thresholds from Section 4.
- Mixed samples (contamination, co-infection, host DNA) are the norm, so the output is a list of matches with abundances, not one label.

Realistic delivery:
- **CLI/HPCC first**, using the branchwater plugin (`fastgather`/`fastmultigather` or the current equivalents) against the T1/T2 indexes. This is mostly integration work.
- **Browser:** sketching a multi-GB FASTQ in wasm is impractical. A *subsample* mode (first ~1–5M reads, streaming, local file) may work and gives a coarse genus/species call. **[est]** Prototype and measure it. It is a stretch goal, not a promise.
- Needs its own benchmark: simulate reads from held-out genomes at 1x/5x/20x with realistic error profiles, and also test real data with known species.

## 6. Goal 3 (PopSketch trees)

Aim: from assemblies, quickly produce a distance matrix, tree and clusters with consistent method choice, and a viewer.

What the pilot showed **[measured]** (15 *A. fumigatus*):
- mash and sourmash gave the same topology (cophenetic Spearman 1.000). skani agreed at 0.88, mainly differing in the placement of eAF1436 and Af293. skani is more resolving for near-identical pairs (A1160/A1163).
- Genomes span ~99.2–99.99% ANI, which is at the edge of what MinHash can resolve.

Constraints to state in the docs:
- Distance-based trees on sketches are **clustering and exploratory visualizations, not phylogenies**. Recombination (common in *A. fumigatus*), accessory content and assembly quality all distort them. Validate against a SNP-based or core-genome tree for a subset (e.g., `parsnp`/`kSNP`/mapping-based) before anyone cites a sketch tree as evolutionary.
- Outbreak or isolate-level resolution (>99.9% ANI) needs SNP calling, not sketches.
- Existing tools to benchmark against: `mashtree`, `skani triangle` + NJ, `sourmash compare` + clustering, `GTDB-Tk`-style ANI clustering. Decide whether PopSketch adds a wrapper, a viewer, or new methods.

From reads: sketch distances between read sets need abundance-filtered sketches and similar coverage. Possible but second-order. Defer until GmG-reads exists, because it reuses the same sketching code.

Proposed default: skani for the final ANI matrix (fast, accurate above 98% ANI) and sourmash FracMinHash for anything that must also run in the browser. Report both methods' agreement as a QC metric.

## 7. MATPredict integration

- Interface: `guess(genome_or_reads) -> {group_id, taxon, confidence, nearest_refs[], flags[]}` as JSON, callable from the CLI and from a library.
- MATPredict uses `group_id` to select the reference MAT locus/gene set. The set of groups MATPredict supports should define a minimum list of clades with high-confidence calls. **Ask: which groups does MATPredict need first?** That prioritizes which clades to benchmark deepest.
- If GmG returns "novel" or low confidence, MATPredict should fall back to a broader reference set instead of silently using the wrong group.

## 8. Milestones with go/no-go gates

| # | Milestone | Gate (decide before moving on) |
|---|---|---|
| 0 | **Measure.** On the 15 genomes plus a broader set (e.g., 100–200 from several genera): sketch at scaled 1000/10000, measure delta-coded and 32-bit-truncated sizes, and species-level accuracy at each. Check genome-length distribution of the full set on HPCC. | Does scaled=10000 still call species correctly? Are the size estimates in Section 3 right (±30%)? |
| 1 | **HPCC build pipeline.** SLURM array sketches all genomes, manifest with taxonomy and provenance, dereplication and clustering. | Reproducible from one config. Tiers T0/T1 produced. |
| 2 | **Benchmark harness.** Leave-one-out and leave-one-species-out, confusion by clade, novelty tests, contaminated and fragmented assemblies. | Accuracy and threshold calibration vs. Section 2 targets. Choose the final scaled/hash width. |
| 3 | **GmG CLI.** Containment/ANI search and taxon aggregation against T0/T1, JSON output. | Meets accuracy criteria. MATPredict can consume it. |
| 4 | **Browser prototype.** wasm sketching + search against T0 on Pages. Check the sourmash Rust crate's wasm support first and don't assume it. | Latency and first-load size within budget. If wasm sketching is too slow, fall back to "upload a precomputed sketch" (sketch locally with the CLI, search in the browser). |
| 5 | **GmG-reads (HPCC).** Branchwater-style search against T1/T2, abundance-filtered reads sketches, simulated and real benchmarks. | Detection limits (min coverage) documented. |
| 6 | **PopSketch trees.** Pipeline + viewer, compare against a SNP tree on a subset. | Concordance documented. Caveats written into the docs. |
| 7 | **(Stretch) Browser reads subsample mode.** | Only if Milestone 5 works and the subsample accuracy is acceptable. |

## 9. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| Browser DB size or hosting limits | GmG can't be client-only | Tiered DB, range requests, object store, "upload sketch" fallback |
| scaled=10000 too coarse for close species complexes | Wrong species-level call within complexes (*Aspergillus* section *Fumigati*, etc.) | Tier T2 refinement by cluster, per-clade thresholds |
| Mislabeled NCBI references | Systematically wrong calls | Top-N voting, ANI-based outlier detection against clustered neighbors, flag disagreements |
| Fragmented or contaminated query assemblies | Low containment, mixed hits | Containment (not Jaccard), contamination flag, minimum-length warnings |
| Novel species not in reference | Forced wrong label | Explicit novelty threshold, report best ANI |
| wasm performance or sourmash wasm support unverified | Browser tool may not be feasible as imagined | Milestone 4 prototype with the fallback above |
| Overlap with existing tools | Duplicate effort | Section 4 prior-art check before Milestone 3 |

## 10. Open questions

1. Exact definition of "species group" and the target clades (see MATPredict needs).
2. Is the ~20k BFD set the full reference, and are the HPCC copies already QC'd and dereplicated?
3. Acceptable first-load size and latency for the browser tool.
4. Who is the user (lab-internal vs. public)? Public use raises hosting, versioning and support costs.
5. Licensing/provenance of shipped sketches (sketches of public NCBI assemblies are generally fine, but record accession and version for each).
6. ~~Does PopSketch live as one repo with GmG?~~ Resolved 2026-10-27: one new repo, GmG. PopSketch (trees/PCA) is a later layer on the same core; revisit splitting it out once it has its own needs.

## 11. Data and artifacts already in hand

- `data/genomes/`: 15 *A. fumigatus* assemblies (accessions in `data/accessions.txt`, candidate table of 416 in `data/afum_candidates.tsv`).
- `data/out/`: mash, sourmash and skani outputs, and `dendro_compare.png` (labelled by strain).
- Conda env `popsketch` (name kept from before the rename): sourmash 4.9.4, mash 2.3, skani 0.4.0, ncbi-datasets-cli, fastani.

## 12. HPCC findings (inspected read-only via `ssh hpcc2`, 2026-10-27)

Base: `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD`

- **`input_clean_genomes/clean/` holds only 2 test genomes** (`GCF_TEST001/002.purge.fasta` + FCS-GX taxonomy TSVs). It is not the full set, and `samples.csv` ASMIDs (e.g. `GCA_022813625.1_ASM2281362v1`) do not match those file names. Confirm which directory is the intended source.
- **The full set appears to be `input/dna/`**: 23,448 symlinks named `{Species}_{Strain}.scaffolds.fa`, pointing into `Fungi_BFD_runs/genome_annotation/*/predict_results/`. Those are funannotate-predict scaffolds (post-cleaning and possibly repeat-masked). Check whether they are soft-masked, because lowercase masking and any renaming can affect sketches. Link name to ASMID needs a join (`{SPECIES}_{STRAIN}` is built in README_BFD.md, ~line 348).
- **`samples.csv`: 23,680 rows**, no duplicate ASMID. Taxonomy is filled in: PHYLUM 23,677, GENUS 23,339, SPECIES 23,680 (7,450 distinct), STRAIN 22,945. 19 phyla, 70 classes, 195 orders, 1,671 genera.
- **Strong skew** (matters for the reference design): 5,576 species have a single genome, 534 species have >=5. Top species: *S. cerevisiae* 1,309; *Candidozyma auris* 952; *F. oxysporum* 814; *P. oryzae* 495; *M. phaseolina* 466; *A. fumigatus* 416. Ascomycota is 19,412 of 23,680. Dereplication will shrink T0 a lot (T0 is at most ~7,450 species reps, and fewer if clustered), and the 5,576 singletons have no within-species replication for leave-one-out testing. Use leave-one-species-out for those.
- **Existing sketching and ANI work to reuse, not redo:**
  - `results/mash_sketch/`: per-genome `.msh` files (12 GB; the file names use ASM-style IDs, not the `input/dna` names). Check mash parameters (`-k`, `-s`) with `mash info` before reuse. It is likely mash-only, with no FracMinHash/containment.
  - `README_ANI.md` and `nextflow/main.nf --pipeline compare_ani`: skani (default), mash, sourmash and fastANI all-vs-all within taxon groups, with a `storeDir` sketch cache (`work/ANI/sketch_cache/` does not currently exist). `results/ANI/skani/` exists. LZ-ANI was evaluated and rejected (see `analysis/ani_method_evaluation/`). The cluster-level and outlier outputs in that pipeline overlap with PopSketch trees (Section 6) and with dereplication (Section 3).
  - `results/query_unknowns_mash/`: an existing "identify unknowns by mash" workflow. Read this first; it overlaps directly with GmG.
- **Actions before Milestone 0:** (1) read `query_unknowns_mash` and `analysis/ani_method_evaluation`; (2) decide the genome source (`input/dna` vs. cleaned set); (3) build the `input/dna` -> ASMID join and check that every row of `samples.csv` resolves to a file; (4) run `mash info` on a few existing sketches.

### 12b. Corrections after checking `Fungi_BFD_runs` (2026-10-27). These supersede the conflicting points in Section 12.

Base: `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs`

- **Genome source = `input_clean_genomes/clean/`** (47,753 files, i.e. two per genome): `{ASMID}.purge.fasta.gz` plus `{ASMID}.purge.fcs_gx-taxonomy.tsv.gz`. About 22k unique assemblies. The `.0`/`.1` suffix variants (1,082 and 187 files) look like split/duplicate outputs, so check them. One earlier listing also showed `.fa.gz` and `.masked.fasta.gz` names, so confirm the set (purged vs. masked) before sketching. Don't sketch soft-masked FASTA unless intended (sourmash is case-insensitive, but be explicit).
- **File names match `samples.csv` ASMIDs here** (`Fungi_BFD_runs/samples.csv`, 23,680 rows). The earlier name-join problem applied to `Fungi_BFD/input/dna` only. Still verify that every row resolves (23,680 rows vs. ~22k assemblies).
- **Existing mash sketches: `results/mash_sketch/`, 22,815 `.msh` files, ~78 KB each (~1.7 GB total).** Pipeline defaults are `k=21, s=10000`, which is consistent with the ~80 KB size (10,000 x 8 B). This is a fixed-size bottom-k sketch, **not** FracMinHash. A bottom-k sketch of 10,000 on a ~30 Mb genome is about the same density as scaled~3000, but sketches from different genome sizes are not on the same hash threshold, so containment of a small query in a large reference is not valid. Good for mash-distance species-level GmG screening at no build cost. Not a substitute for scaled sketches (T0/T1/T2). Verify with `mash info` (mash is not on the login-node PATH; use the conda env or the pipeline's container).
- **Existing skani results: `results/ANI/skani/{GENUS,SPECIES}/<taxon>/`**, with per-taxon `*_ANI_report.txt`, `*_genome_names.tsv` and `batches/` (935 genus dirs, 1,871 species dirs; run Aug-Oct 2025). Within-taxon only (skani `medium` preset, min_af=15), so there is no cross-genus matrix. Use these as ground-truth species clusters and thresholds for GmG calibration and for dereplication, instead of recomputing.
- **Pipeline defaults for sourmash: `k=21, scaled=1000`** (`nextflow.config`, `nextflow/modules/ani/sketch/`). Check whether sourmash sketches were ever materialized (the `work/ANI/sketch_cache` path did not exist under `Fungi_BFD`; look under `Fungi_BFD_runs/work`). If not, producing them is Milestone 1.
- **Revised Milestone 0/1 scope:** (a) reuse the existing mash sketches for a quick species-level baseline to see how far mash-only GmG gets (leave-one-out using `samples.csv` SPECIES labels and the skani species clusters); (b) build FracMinHash scaled=1000 + 10000 from `input_clean_genomes/clean` in one pass (~7 core-hours **[est]**), since only that supports containment, `gather` and the tiered browser DB.

### 12c. Genome source corrected again (2026-10-27). Supersedes the genome-source bullet in 12b.

- **Use `Fungi_BFD_runs/input_clean_genomes/` directly, not `.../clean/`.** It holds `{ASMID}.fa.gz` (cleaned, unmasked) and `{ASMID}.masked.fasta.gz` for each genome: 47,746 files, so about 23.9k genomes by file count. That is slightly more than the 23,680 rows in `samples.csv`, so the join needs a check in both directions (files with no row, rows with no file). `clean/` is a separate subdirectory that I looked at by mistake.
- **The `.0/.1/.2` "variants" I flagged in 12b are not variants.** They are the assembly-version digits that end some ASMIDs (e.g. `..._1.0.fa.gz`) and were an artifact of my filename parsing. No action needed.
- **Sketch the unmasked `.fa.gz`.** The `.masked.fasta.gz` copies are for repeat-masked annotation, and masking would remove repeat k-mers that are useful for identification. Whether masking is soft or hard is unchecked, so don't use them unless there is a reason.
- The `.fcs_gx-taxonomy.tsv.gz` files are in `clean/` only; if FCS-GX contamination calls are wanted as a QC flag, look there.

### 12d. Containers and tool versions in the existing BFD ANI pipeline (2026-10-27)

Source: `Fungi_BFD/nextflow/conf/profile_ANI.config`, `modules/ani/sketch/*/main.nf`. Cache dir: `/bigdata/stajichlab/shared/singularity_cache`.

| Tool | Version | Container (pre-built sandbox) |
|---|---|---|
| mash | 2.3 | `mash_2.3--hb105d93_10.sandbox` (BioContainers) |
| sourmash | 4.9.4 | `sourmash_4.9.4--hdfd78af_0.sandbox` |
| skani | 0.3.2 | `skani_0.3.2--h79ce301_0.sandbox` (0.2.2 `.img` also cached) |
| fastANI | 1.34 | `fastani_1.34--hb66fcc3_5.sif` |
| report/python | 3.12 | `python_3.12.sandbox` |

- Pipeline defaults: mash `-k 21 -s 10000`, sourmash `k=21,scaled=1000`, skani preset `medium` in `nextflow.config` but `fast` (`-c 200`) in `profile_ANI.config`, so confirm which preset produced `results/ANI/skani/`. The sketch and compare steps must use the same preset.
- `MASH_SKETCH` is `mash sketch -k 21 -s 10000 -o {genome.baseName} {genome}`, with `storeDir`, 2 CPU / 4 GB, one genome per task. This confirms the existing `.msh` files are bottom-k at the pipeline defaults (the sketch files are ~80 KB, consistent with 10,000 x 8 B). I have not run `mash info` on them (see below).
- **Use the same container versions in PopSketch for the HPCC build** (sourmash 4.9.4, skani 0.3.2, mash 2.3) so that new sketches are comparable with the existing ones. My local test env uses skani 0.4.0 and the same sourmash/mash versions, so re-check skani agreement on the 15-genome set with 0.3.2 before mixing results.
- Containers are pre-built **sandboxes**, not `.sif` or `docker://` URIs, on purpose: squashfuse timeouts and full node `/tmp` broke high-concurrency runs on 2026-08-26/27 (per comments in the config). Follow the same convention for any new SLURM array. Do not pull `docker://` images at job time.
- `apptainer`/`singularity` is not available on the login node (no `module` either), so `mash info` on the existing sketches needs a compute-node allocation. **Pending**: run `apptainer exec -B /bigdata <mash sandbox> mash info -H <one .msh>` inside a short interactive job to verify k and s.

### 12e. Verified (2026-10-27): resolves the "Pending" item in 12d

- Running containers on HPCC: the login node has no `apptainer` binary even after `module load apptainer` (loads 1.5.3 but only on compute nodes). Use `srun -p short ... bash script.sh` with `module load apptainer` inside the script. Build jobs should do the same. Nested-quote `srun bash -c "..."` through ssh mangled the variables, so use a script file.
- `mash info` on `results/mash_sketch/GCA_000150145.1_ASM15014v1.fa.msh`: **k=21, 10,000 min-hashes, MurmurHash3_x64_128 seed 42, canonical, 64-bit hashes**, sketched from `GCA_000150145.1_ASM15014v1.fa.gz` (29,174,577 bp, 32 seqs). So the existing sketches were made from the unmasked cleaned `.fa.gz`, which matches the recommendation in 12c.
- Container versions confirmed inside the sandboxes: sourmash 4.9.4, skani 0.3.2, mash 2.3.
- Still open: whether the skani results used `fast` or `medium`.

### 12f. skani preset provenance (checked 2026-10-27)

- **Cannot be verified from what is on disk.** The `results/ANI/skani/GENUS/*` outputs are dated 2026-08-27 and the `SPECIES/*` ones 2026-10-07. The ANI reports record the cluster/outlier thresholds but not the skani preset. The `.nextflow.log*` files in `Fungi_BFD_runs` only start 2026-09-05, and no `work/` task dirs remain from the Aug 27 run.
- `profile_ANI.config` (`skani_preset = 'fast'`) was last modified **2026-08-28 10:41**, the day after the GENUS run, so it may not have been in effect then. The base `nextflow.config` and `README_ANI.md` say `medium`. The SPECIES run (Oct) probably used the current `-profile ani` settings (`fast`), but that is an inference.
- **Decision:** treat the existing skani ANI as "preset unknown, `fast` or `medium`". Use it for species-cluster ground truth (95% threshold) and as a sanity check, not as a precision benchmark. For anything PopSketch publishes, re-run skani with an explicit, recorded preset on the subset that matters. Record `skani_preset`, `-c`, `min_af` and container version in every result manifest.

### 12g. skani preset unified (2026-10-27)

- `Fungi_BFD/nextflow/nextflow.config` line 63 changed `skani_preset` from `'medium'` to `'fast'` to match `conf/profile_ANI.config` (`min_af = 15` was already identical). Backup: `nextflow.config.bak_20261027` next to it. Nothing re-run.
- Still `medium`, left untouched: `conf/profile_funannotate.config:559`, `k8/profile_compare_ani_k8s.config:104`, and the default documented in `README_ANI.md:106`. Revisit if those pipelines should also use `fast`.
- Provenance caveat from 12f still applies to existing results (Aug GENUS run preset unknown).

## 13. Decisions and new findings (2026-10-27). These supersede the open questions in Section 10.

### 13.1 Decisions from the lab
1. **Species group:** ship both an NCBI-taxonomy label and an ANI-cluster ID.
2. **MATPredict needs phylum first:** Ascomycota, Basidiomycota, Mucoromycota. Whether a specific NCBI taxid is also needed is addressed in 13.2.
3. **Reference/representatives from empirical data**, using the BFD ANI-clustering approach (13.3), with NCBI representative genomes only as a cross-check. NCBI "representative" selection hides complexity (mislabels, multiple strains per species).
4. **Size:** target ~100 MB for the default browser payload. Up to ~2 GB is acceptable for an optional full download, but prefer representatives for speed.
5. **Scope now:** a new tool. Build the Python/Rust tool interface first (library + CLI). Deployment, command line and WASM are decided afterwards.

### 13.2 MATPredict contract (read from `~/projects/MATPredict`, v0.7.0)
- `matpredict detect --genome G [--proteins P] [--taxid TAXID] [--phylum {Ascomycota,Basidiomycota,Mucoromycota}]`. **`--taxid` is optional.** The taxid selects the curated families by lineage (and the genetic code). Without it, MATPredict uses `--phylum` (`explicit_phylum`), falls back to a phylum guess (`phylum_fallback`), or reports `not_searched`. Other routing modes: `lineage`, `exhaustive`.
- **So the minimum GmG output for MATPredict is a phylum call, and the more useful output is an NCBI taxid** (nearest-reference `NCBI_TAXONID` from `samples.csv`) so MATPredict can route to `lineage` mode. Both go in the JSON result with a confidence flag. If the confidence is low, GmG should return only the phylum, or none. MATPredict falls back to `not_searched`/`exhaustive` instead of searching the wrong families. **Open:** does MATPredict need the taxid to be species-level or is order/family enough? (`db/<Phylum>/<Order>/` routing suggests order/family is the working level; check `taxonomic_scope` in the family `order.yml` files.)
- MATPredict already does **reads-based work**: `datasets/reads/sra_survey/pilot/` uses a sourmash **k=21 FracMinHash sketch with abundances of a read subsample, `SKETCH_SCALED=500` (~100 kB per run)** (`genotype_mat_one_accession.sh`, `clone_clusters.py`), calling the HPCC sourmash by path (`/opt/linux/rocky/8.x/x86_64/pkgs/sourmash/4.9/bin/sourmash`). **Implication: make GmG-reads compatible.** Sourmash can downsample scaled=500 to 1000 but not the reverse, so build the references at scaled <= 1000 and accept query sketches at scaled <= 1000. Reuse the subsample convention instead of inventing a new one.

### 13.3 BFD representative/cluster logic (read from `Fungi_BFD/nextflow/bin/`)
- `pick_representative_strain.py` / `species_reuse_clusters.py`: per species with >=2 strains, the representative is the strain with the **highest BUSCO genome completeness, then highest scaffold N50, then alphabetical** as the final tiebreak. Other strains are `reuse_eligible` if skani ANI to the representative is **>= 99.0%**, and missing ANI is treated as ineligible (fail-closed). Inputs: `tables/busco_genome.parquet`, `tables/asm_stats.parquet` (keyed by ASMID), and the merged ANI table `results/ANI/skani/SPECIES/all_pairs_merged.tsv` (or `ani.duckdb`).
- Thresholds in `conf/profile_ANI.config`: `ani_cluster_threshold = 95.0`, `ani_outlier_threshold = 90.0` (used by `scripts/find_ani_label_mismatches.py` to audit cluster-vs-label disagreements). **`find_ani_label_mismatches.py` is directly useful for GmG:** it finds the mislabeled references that would otherwise poison species calls. Use its output to drop or flag those references.
- **Proposed T0/T1 derivation:** T1 = representatives at the **99%** BFD rule (near-clonal strains collapsed). T0 = one representative per **95% ANI cluster** within each species, chosen with the same BUSCO > N50 > name rule, so the choice is auditable and shares code with BFD. The 99% value in BFD was chosen for ab-initio parameter reuse, not for identification, so it may be too strict or too loose for T0. Test it on the benchmark. Large species (*S. cerevisiae* 1,309, *C. auris* 952, *F. oxysporum* 814) will collapse to a handful of reps, while the 5,576 singleton species stay as they are.
- Note: skani `fast` vs `medium` still unverified for the existing matrices (see 12f). Use the 95% cluster membership as a start and re-check borderline pairs.

### 13.4 Existing GmG precursor: `analysis/desert_microbes_mash_id`
- Already identifies 14 unlabeled desert-isolate assemblies against the BFD genomes: **`mash dist` against 23,867 sketched genomes (k=21, s=10000), ~49 min on 16 CPUs**, results joined to `samples.csv`. 10 of 23,877 files had zero FASTA records and had to be excluded (`mash sketch -l` aborts the whole run on one empty file, so validate first). The `.fa.gz` count is therefore 23,877, which resolves my earlier "~23.9k" estimate.
- Its own findings set the expectations for GmG: a near-species match (98.2% ANI, 5,265/10,000 shared hashes, *Aureobasidium melanogenum*) works well, but mash ANI is only a rough proxy and **gets noisy below ~80%**, where shared-hash counts fall to single digits. Several isolates resolved only to family or class or to Basidiomycota with very weak signal. **GmG must report a confidence tier, and "phylum only" and "unknown" are valid answers.**
- `results/query_unknowns_mash/` (`candidate_lineages.tsv`, `lineage_triage.html`, `query_list.tsv`) is a related workflow. Not yet read.
- **What GmG adds over this:** one reusable tool, not a one-off analysis. It adds containment/FracMinHash (works on incomplete assemblies), a precomputed small index, calibrated thresholds validated by leave-one-out, MATPredict-ready output, and (later) reads and browser delivery.

### 13.5 Updated milestone sequence
| # | Step | Notes |
|---|---|---|
| 0a | Reuse existing mash sketches (22,815 `.msh`) to do a leave-one-out phylum/order/species baseline against `samples.csv` | no new sketching; compare against 13.4's noise limit |
| 0b | Build the cluster/representative table: `all_pairs_merged.tsv` + `busco_genome.parquet` + `asm_stats.parquet` + label-mismatch audit | check output of `find_ani_label_mismatches.py` |
| 1 | Define the result schema (JSON) incl. phylum, taxid, ANI/containment, confidence tier, routing hint for MATPredict | the contract with MATPredict, reviewed with it |
| 2 | Build the core library + CLI (Rust core with Python bindings, see 13.6), reading a reference index in the agreed format | no wasm yet |
| 3 | Benchmark harness (leave-one-out, leave-one-species-out, novel-species, contaminated or fragmented assemblies) with phylum-level accuracy as the first gate for MATPredict | |
| 4 | Scaled sketches (k=21, scaled <= 1000) from `input_clean_genomes/*.fa.gz` for the representative set first (not all 23.9k), then the rest if the benchmark needs them | validate FASTA records before sketching; use the BFD container sandboxes via SLURM |
| 5+ | Reads (GmG-reads), PopSketch trees, then WASM/deployment (decided after the CLI works) | |

### 13.6 Interface-first proposal (to be confirmed before coding)
- **Core in Rust, Python bindings (PyO3), CLI on top; WASM later from the same core.** Reason: the lab already has Python tooling (BFD, MATPredict), HPCC needs a fast native path, and a Rust core is the only route to wasm without a rewrite. **Verify before committing:** that the `sourmash` Rust crate (or a thin FracMinHash implementation compatible with sourmash 4.9.4 sketches) builds for the targets we need. I have not checked this, and local `rustc`/`cargo` is available but `wasm-pack` is not installed.
- **Alternative if Rust is too heavy at first:** a pure Python CLI on the sourmash Python API for the reference build and HPCC use, with the Rust/WASM path developed in parallel for the browser. Decide after the result schema (step 1) so both can share it.
- Proposed commands: `gmg build` (HPCC: sketches + clustering -> index), `gmg guess GENOME` (search, JSON out), `gmg reads FASTQ...` (later), `gmg info INDEX`.

## 14. Rust core: decision, evidence, first code (2026-10-27)

### 14.1 Decisions (lab)
- **Rust core first** (Section 13.6 resolved).
- **Taxon level:** MATPredict needs order only. **GmG is a general tool, so aim for genus and species and measure whether that works** (benchmark per rank: phylum / order / genus / species, with a confidence tier per rank and "stop at the rank we can support").

### 14.2 Feasibility evidence (scratch tests, not kept)
- The `sourmash` Rust crate (0.23.0, 2026-06-07) reproduces Python sourmash 4.9.4 exactly on A1163 vs Af293 at k=21, scaled=1000: same hash counts (28,608 / 28,276), same intersection (26,329) and union (30,555), same first hashes. Sketching ~0.4 s per genome natively vs ~1.2 s with Python sourmash (A1163 + Af293 in 0.8 s total, one thread).
- **The `sourmash` crate is a poor base for the wasm build.** It compiles for `wasm32-unknown-unknown`, but it builds its own `cdylib` C-FFI, which fails to link because it drags in `zstd-sys` C code, and it has a heavy dependency tree (zstd, lzma, flate2, csv, nalgebra, memmap2, rayon optional). **Decision: write a small, dependency-light FracMinHash in `gmg-core` and keep the `sourmash` crate as a dev-only oracle.**
- The rustup `stable` toolchain's `rust-lld` aborts here with a missing `libLLVM.dylib`, and the Homebrew `rustc` (first on PATH) has no wasm std. Building wasm needs `PATH=$HOME/.rustup/toolchains/stable-aarch64-apple-darwin/bin:$PATH`. The lld problem means a full `cdylib` wasm link has **not** been tested yet (only `rlib` builds). Fix or work around (e.g. `wasm-ld` from Homebrew LLVM, or `wasm-pack`/`trunk` once installed) when the wasm step starts.

### 14.3 What is built
Workspace at the project root: `Cargo.toml`, `crates/gmg-core/`.
- `gmg-core` (no I/O, deps: `murmurhash3` only): `FracSketch` (k, scaled, seed; case-insensitive; skips non-ACGT; canonical k-mers; sourmash `max_hash`), `containment()`, `jaccard()`.
- Tests (all passing): 3 unit tests, plus `tests/oracle.rs`, which checks **hash-for-hash equality with the `sourmash` crate on the A1163 and Af293 pilot genomes** (skips if `data/genomes/` is absent).
- `cargo build --release --target wasm32-unknown-unknown -p gmg-core` succeeds (rlib).
- `.gitignore` added for `target/` and the downloaded data. The project is still not a git repo.

### 14.4 Next code steps (proposed order)
1. **Result schema** (JSON, per-rank calls + confidence + MATPredict routing hint). Needs a quick check of MATPredict's order-level `taxonomic_scope` field names so the hint can be consumed directly.
2. **Index format v0** in `gmg-core`: header, manifest, sorted hash blocks, optionally delta-coded, with an offset table for range requests. Round-trip tests and size measurement on the 15 genomes at scaled=1000/10000.
3. **`gmg-cli`** crate: `gmg sketch` (FASTA/gz via `needletail`, native only), `gmg guess`, `gmg info`.
4. **Python bindings** (PyO3/maturin) once the API settles, so MATPredict and BFD scripts can call the core.
5. **Benchmark harness** against the existing mash sketches (Milestone 0a) to see how far genus and species calls go before the scaled index is built.
