# Web app

`web/` is a static page that runs GmG in the browser: the same Rust code as the `gmg` command line, compiled to WebAssembly
(`crates/gmg-wasm`). Files are sketched on your computer and never uploaded. Results (records, `results.tsv`, `hits.tsv`
and the report) are identical to `gmg run` for the same inputs and index (floating-point values agree to the last digit or
two); `crates/gmg-wasm/tests/parity.mjs` checks this on 21 inputs (assemblies, mutants, reads, empty and invalid files) in CI.

* **Index:** pick a bundled index (the representatives index is published with the app) or any `.gmgi` file, for example the
  full index from the [index releases](https://github.com/stajichlab/GmG/releases). The index is cached in the browser
  (Cache API), so the next visit works offline.
* **Inputs:** drop files or a folder (`.fa .fna .fasta .fas .ffn .fq .fastq`, plain or `.gz`; gzip is decompressed in the
  browser). FASTQ: only the first N reads are read (default 500,000). The work runs in a Web Worker, with progress and Cancel.
* **Outputs:** a sortable, filterable table; downloads of `results.tsv`, `hits.tsv` and `report.html`; "Print / save as PDF"
  opens the report in its own tab for the browser's print dialog.
* **Not in the browser:** merging R1/R2 and lane files into one sample (each file is a row), sample sheets, `--resume`.

## Build and run locally

```
pixi run wasm                                   # builds crates/gmg-wasm/pkg-web
scripts/build_site.sh _site bfd-reps-2026.10.1-s20000.gmgi
python3 -m http.server -d _site 8000            # the worker needs http(s); file:// does not work
```

## Publishing

`.github/workflows/pages.yml` builds the site on every push to `main` that touches the web app or the Rust code, bundles the
representatives index of the newest `index-*` release and deploys to GitHub Pages (repository setting: Pages source =
GitHub Actions) and the repository variable `PAGES_ENABLED=true`; until then the workflow only runs by hand. The full index (200 MB) is not bundled: load it from a file.

Memory: the browser holds the index once (about its file size plus the sketch of the current file): ~180 MB for the
representatives index, ~460 MB for the full index (`reports/performance_2026-10-28.md`).

Tested by hand in headless Chromium (load a bundled index, classify genomes and a FASTQ.gz, table, downloads, report, no console
errors). Firefox and Safari have not been tried; the page uses ES modules in a module worker, `DecompressionStream`, the
Cache API and `<input webkitdirectory>`, which current versions of all three support.
