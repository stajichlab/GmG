# Changelog

## Unreleased

- **Marker mode** (`docs/marker-mode.md`): identify a *Fusarium* isolate from marker sequences (tef1, rpb1, rpb2, CaM, tub2, ITS) against a FusariumID strain index built at scaled 1 (`scripts/fusariumid_markers.py prepare`, then `gmg sketch -s 1` and `gmg build`). New built-in threshold set `--thresholds markers` (species ANI 0.995), calibrated on a held-out 20% of FusariumID strains (`bench/F2_marker_mode/`). The index itself is not distributed (FusariumID data).
- **Species complexes (schema 10).** Results carry `species_complex` (JSON `guess.species_complex` with `complex`, `nearest_species`, `collapsed`; new TSV column before `schema_version`): the species vote's complex from an embedded table of 806 fusarioid species and informal taxa built from the FusariumID database (`scripts/make_species_complexes.py`), reported when the genus is called and the top hit is within ANI 0.98. For the *F. oxysporum* species complex the species call is withheld (`withheld` code `species_complex`) and the complex is reported instead: its species are 99.6-99.99% genome ANI apart and 38% of FOSC genomes have another FOSC species above the species cutoff (`bench/F1_fusariumid/`). Other species calls are unchanged. `expected_match` also accepts the complex name.
- **popsketch** (new crate and binary, `docs/popsketch.md`): `popsketch run reads/ -o out/ --labels labels.tsv` audits a set of read files and shows coarse population structure. The audit flags, with the reason and numbers for each, samples that are low in coverage, far from the cohort in estimated size, carry an abundance shoulder (a second organism at several-fold coverage), or share almost nothing with the other samples; with `--index` it also names the organism and finds mixtures. Samples that pass get a hash-presence PCA, clusters (cut only where the jump in merge height is at least 3 measured noise floors), a neighbour-joining tree with branches below the noise floor collapsed, an HTML report, and a comparison with a multi-column label table (adjusted Rand index; samples whose ANI matches another label better than their own). `--flagged exclude|mark|decontaminate`, `--max-samples` (default 1500). Light contamination (about 1 to 7% of reads) is not visible without `--index`.
- **popsketch**, further options: `--compare-tree` and `--compare-pcs` compare a run with a SNP tree and PCA of the same samples (distance correlation, whether each cluster is a clade, principal angles); DAPC on the primary label column (`--dapc-pcs`, `dapc.tsv`); `--report-format html|md|pdf|all` (Markdown with standalone SVG figures, and a PDF); `--target-cov` stops reading a sample once its k-mer coverage reaches the target; exact all-pairs distances by bitsets (identical result, about 4x faster at 300 samples; 1,500 samples take about a minute after sketching, 3,000 about 4.5 minutes); `popsketch demo` writes a small simulated example and `docs/popsketch-tutorial.md` walks through it. The abundance-shoulder flag is measured after removing the part that tracks heterozygosity, and is a note and not an exclusion when more than 2% of a cohort has one.
- New library crate `gmg-input`: finding FASTA/FASTQ files, naming samples, merging R1/R2 and lanes, and reading sample sheets, moved out of `gmg run` without a change of behavior so `popsketch` shares the rules.
- Read sets: the abundance cutoff (`auto_min_abund`) and k-mer coverage no longer take a spurious valley from an error-dominated histogram. A valley is accepted only if the hashes at or above it are at least half of the hashes seen twice or more. Before, a read set with essentially no coverage peak (sub-1x, or a low-coverage plateau with no rise after it) could get a cutoff of 15 to 29 and a "coverage" of 23 to 420x, which left 0.1 to 3 Mb of hashes. On the 559 Coccidioides and 334 *B. dendrobatidis* read sets used to check this, exactly four samples change (to cutoff 2 and a coverage of 0.7 to 1.3x) and none of the others do. Calls on read sets with a real peak are unchanged.

## v0.3.2 (2026-11-04)

Calls and schema 9 are unchanged from v0.3.1; the change is which files are read as sequence. The Python wheels (PyPI `guess-my-genome`) carry it.

- FASTA files ending in `.fsa` (the GenBank submission name), `.seq` or `.dna` are found and classified like `.fa`/`.fna` (command line, folders, and the browser app); the sample id drops the extension. Before, `gmg run` answered "no FASTA/FASTQ files found" for them. The format is read from the content, so a GenBank-format `.seq` file is reported as `not_sequence_file`, not misread.

## v0.3.1 (2026-11-03)

No change to the program or its results (same calls and schema 9 as v0.3.0); the Python bindings are published to PyPI for the first time.

- PyPI: `pip install guess-my-genome` (import name `gmg`). A release workflow (`.github/workflows/pypi.yml`) builds wheels for Linux (x86_64, aarch64) and Apple Silicon and publishes them with trusted publishing when a version tag is pushed; richer package metadata. Wheels only: there is no source distribution, because the build reads `data/thresholds/` from outside the crate; other platforms build from the repository with maturin. Intel Macs and Windows may follow.

## v0.3.0 (2026-11-02)

Built to run with `index-2026.10.5` (no new index; thresholds unchanged, calls are identical to v0.2.1).

**For tools that route by phylum (MATPredict)**
- Results schema 9: every record (`results.jsonl`, Python `Index.classify`, WASM) has a `routing_hint`: `usable`, `reason`, `phylum`, `phylum_taxid`, `deepest_rank`, `deepest_taxid`. `usable` is false when there is no phylum call, the input is flagged `not_fungal`, or the call is `weak_phylum`; the call is still reported so a consumer can apply its own rule. The TSV columns are unchanged (only `schema_version` reads 9). See `docs/results-schema.md`. Measured on MATPredict's regression panel: no wrong usable call when relatives are in the reference; with the genus hidden, 4 of 40 usable calls were wrong, all into a phylum MATPredict has no families for.
- Python bindings: the PyPI distribution is named `guess-my-genome` (`gmg` on PyPI is an unrelated project); `import gmg` is unchanged. Not published to PyPI yet; build it with maturin or use a release binary.

**Browser app and report**
- The footer shows the GmG version (read from the wasm package).
- Report: in-page links work inside the web app's report frame.

## v0.2.1 (2026-10-31)

- Browser app: GoatCounter page views (anonymous, no cookies; disclosed in the footer) and a footer line "N classifications started so far" fed by a count of presses of Classify (no file information is sent).
- Browser app: loading the index shows a thicker progress bar with MB, percent and time left, then a "building the index" step; the bar was wrong (past 100%) when the server compressed the file, and Load is disabled while loading.
- Browser app: drop or choose a `.zip`, `.tar` or `.tar.gz` (also `.tgz`) and the FASTA/FASTQ files inside are listed and classified like any others (unpacked in memory, up to 3 GB of sequence files per archive).
- Report: the genus and species threshold labels in the ANI histogram no longer overlap.
- Browser app: the full index is offered next to the representatives index.
- Release workflow: the `.dockerbuild` record is no longer attached to the GitHub release.
- The command line reads `.zst`, `.bz2` and `.xz` input as well as `.gz` (folders, `gmg guess`, `gmg sketch`); results are identical. The browser app still reads only `.gz` (the browser decompresses gzip natively).

## v0.2.0 (2026-10-31)

First tagged release. Indexes are released separately (`index-YYYY.MM[.N]`); this version is built to run with `index-2026.10.5`.

**Calls and flags**
- Genus calls need 0.90 vote consistency (was 0.50); family ANI stays 0.7715 (a stricter family margin was tried and reverted, #58).
- QC and call flags: `not_fungal`, `distant_non_fungal_hit`, `weak_phylum`, `mixture_suspected` (flag only, never changes the call, #31), `high_error_reads` (long or noisy reads, #44).
- Results schema 8 (`explained_fraction`, `second_taxon`, `second_explained_fraction`, `qc.coverage_est`).

**Command line and reports**
- `gmg run` folder in, `results.tsv` + `hits.tsv` + `run.json` + HTML/PDF/Markdown report out; reads (DNA, RNA, long) with `--min-abund auto`.
- Compact summary strip at the top of the HTML report.

**Browser app** (gmg.fungalgenomes.org)
- New "Still water" page; the index dropdown names the release; a list of the files to classify (name, size, SHA-256) with a remove button; fetch any NCBI assembly by accession.

**Validation** (reports/validation/)
- Held-out and truth-set precision, independent datasets, stress tests, reads (SRA, mixtures, long reads), mixture detector, comparison with skani, Mash and sourmash; limits in docs/limitations.md.

**Fixes**
- `--version` showed a stale commit after `git pull` (build.rs now watches the branch ref); the HPCC deploy script refuses a binary whose version does not carry the checked-out commit.
