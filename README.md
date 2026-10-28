# GmG: Guess My Genome

Given a fungal genome (or reads), guess where it belongs: phylum, class, order, family, genus, species, plus an
ANI-cluster, each with a confidence. Built for MATPredict routing, meant as a general tool.

Status: working CLI and core (`gmg`), a WASM build, and a 23,669-genome BFD reference built on the HPCC.
Results: [reports/overnight_2026-10-27.md](reports/overnight_2026-10-27.md). History: [LOG.md](LOG.md). Design: [plan.md](plan.md).

## Develop (no Homebrew or rustup needed)

Everything runs in a [pixi](https://pixi.sh) environment (Rust, wasm tooling, Node, Python, optional bio tools):

```
pixi run test        # cargo test (the sourmash-crate oracle test needs data/genomes, see data/accessions.txt; skipped if absent)
pixi run build       # target/release/gmg
pixi run wasm-test   # builds gmg-wasm (nodejs + web bindings in crates/gmg-wasm/pkg-*) and checks it against the CLI
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
| `crates/gmg-wasm` | JavaScript API (`Index`: load, sketch, guess) |
| `scripts/` | analysis and audit (Python), `scripts/hpcc/` SLURM jobs |
| `reports/`, `data/` | results, calibrated thresholds, cluster tables, label audit |

## HPCC

A checkout lives in `~/projects/GmG` (`gh repo clone stajichlab/GmG`); work is committed here, pushed, then `git pull` there.
Build with `pixi run build` in that checkout (login node, short compile) and install with `scripts/hpcc/deploy_gmg.sh`
(versioned binary + symlink swap; never overwrite a running binary). Data and results are under
`/bigdata/stajichlab/shared/projects/GmG_runs/`.
