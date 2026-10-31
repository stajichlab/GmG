# Changelog

## Unreleased

- Browser app: GoatCounter page views (anonymous, no cookies; disclosed in the footer) and a footer line "N classifications started so far" fed by a count of presses of Classify (no file information is sent).

- Browser app: loading the index shows a thicker progress bar with MB, percent and time left, then a "building the index" step; the bar was wrong (past 100%) when the server compressed the file, and Load is disabled while loading.

- Browser app: drop or choose a `.zip`, `.tar` or `.tar.gz` (also `.tgz`) and the FASTA/FASTQ files inside are listed and classified like any others (unpacked in memory, up to 3 GB of sequence files per archive).
- Report: the genus and species threshold labels in the ANI histogram no longer overlap.
- Browser app: the full index is offered next to the representatives index.

- The command line reads `.zst`, `.bz2` and `.xz` input as well as `.gz` (folders, `gmg guess`, `gmg sketch`); results are identical. The browser app still reads only `.gz` (the browser decompresses gzip natively).

## v0.2.0 (unreleased)

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
