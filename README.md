# GmG: Guess My Genome

[![CI](https://github.com/stajichlab/GmG/actions/workflows/ci.yml/badge.svg)](https://github.com/stajichlab/GmG/actions/workflows/ci.yml)

Given a fungal genome (or reads), guess where it belongs: phylum, class, order, family, genus, species, plus an
ANI-cluster, each with a confidence. Built for MATPredict routing, meant as a general tool.

Status: working CLI and core (`gmg`), a WASM build, and a 23,669-genome BFD reference built on the HPCC.
Results: [reports/overnight_2026-10-27.md](reports/overnight_2026-10-27.md). History: [LOG.md](LOG.md). Design: [plan.md](plan.md).

## Quick start

```
gmg run my_genomes/ --index bfd_s10000.v1.gmgi -o results/     # a folder (or files) of FASTA / FASTQ
column -t -s$'\t' results/results.tsv | less -S                 # one row per input; hits in results/hits.tsv
```

`results/report.html` is the report (open it in a browser); `results.tsv` has one row per input with the calls from phylum to species, confidence, QC flags and status; `run.json` records how the
run was made. Columns are documented in [docs/results-schema.md](docs/results-schema.md). Other commands: `gmg sketch / build / guess / info / convert / cluster / bench / completions / man`.
In the browser: [docs/web.md](docs/web.md). Example genomes to try: [examples/](examples/README.md). From Python: [crates/gmg-py](crates/gmg-py/README.md). Step by step: [docs/tutorial.md](docs/tutorial.md); what the results do and do not mean: [docs/limitations.md](docs/limitations.md); every command has examples in `gmg <command> --help`.

## Develop (no Homebrew or rustup needed)

Everything runs in a [pixi](https://pixi.sh) environment (Rust, wasm tooling, Node, Python, optional bio tools):

```
pixi run test        # cargo test (the sourmash-crate oracle test needs data/genomes, see data/accessions.txt; skipped if absent)
pixi run build       # target/release/gmg
pixi run wasm-test   # builds gmg-wasm (nodejs + web bindings in crates/gmg-wasm/pkg-*) and checks it against the CLI
pixi run pytest      # Python tests for the analysis scripts
pixi run fmt         # rustfmt check
pixi run clippy      # lints, warnings are errors
pixi run -e bio ...  # adds sourmash 4.9.4, mash 2.3, skani 0.3.2, ncbi-datasets (same versions as the HPCC containers)
```

`pixi.lock` pins everything for macOS arm64 and linux-64.

## Layout

| path | what |
|---|---|
| `crates/gmg-core` | FracMinHash sketching (sourmash-compatible), index, search, taxon calls, clustering helpers; no I/O |
| `crates/gmg-cli` | `gmg sketch / build / guess / info / convert / cluster / bench` |
| `crates/gmg-report` | report engine: standalone HTML report with inline SVG charts (no I/O; also runs in WASM) |
| `crates/gmg-wasm` | JavaScript API (`Index`: load, sketch, guess; `render_report`) |
| `scripts/` | analysis and audit (Python), `scripts/hpcc/` SLURM jobs |
| `reports/`, `data/` | results, calibrated thresholds, cluster tables, label audit |

## HPCC

A checkout lives in `~/projects/GmG` (`gh repo clone stajichlab/GmG`); work is committed here, pushed, then `git pull` there.
Build with `pixi run build` in that checkout (login node, short compile) and install with `scripts/hpcc/deploy_gmg.sh`
(versioned binary + symlink swap; never overwrite a running binary). Data and results are under
`/bigdata/stajichlab/shared/projects/GmG_runs/`.

## Releases and container

Tagging `vX.Y.Z` runs `.github/workflows/release.yml`: `gmg` binaries (linux-64, osx-arm64), the wasm package and an OCI image,
`ghcr.io/stajichlab/gmg:X.Y.Z` (and `:latest`). Every merge to `main` also publishes `:edge`. Image:

```
docker run --rm -v "$PWD:/data" ghcr.io/stajichlab/gmg:edge guess --index /data/bfd.gmgi /data/genome.fa
apptainer pull gmg.sif docker://ghcr.io/stajichlab/gmg:edge      # on the HPCC (module load apptainer, on a compute node)
```
