# Changelog

## Unreleased

- PyPI: a release workflow (`.github/workflows/pypi.yml`) builds `guess-my-genome` wheels for Linux (x86_64, aarch64) and Apple Silicon and publishes them with trusted publishing on a version tag; richer package metadata. Wheels only (no source distribution; the build reads `data/thresholds/` from outside the crate).

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
