# GmG: Guess My Genome

[![PyPI](https://img.shields.io/pypi/v/guess-my-genome?label=PyPI)](https://pypi.org/project/guess-my-genome/)
[![Python](https://img.shields.io/pypi/pyversions/guess-my-genome)](https://pypi.org/project/guess-my-genome/)
[![Release](https://img.shields.io/github/v/release/stajichlab/GmG?label=release)](https://github.com/stajichlab/GmG/releases/latest)
[![CI](https://github.com/stajichlab/GmG/actions/workflows/ci.yml/badge.svg)](https://github.com/stajichlab/GmG/actions/workflows/ci.yml)
[![Python bindings](https://github.com/stajichlab/GmG/actions/workflows/python.yml/badge.svg)](https://github.com/stajichlab/GmG/actions/workflows/python.yml)
[![Golden panel](https://github.com/stajichlab/GmG/actions/workflows/golden.yml/badge.svg)](https://github.com/stajichlab/GmG/actions/workflows/golden.yml)
[![Fuzz](https://github.com/stajichlab/GmG/actions/workflows/fuzz.yml/badge.svg)](https://github.com/stajichlab/GmG/actions/workflows/fuzz.yml)
[![Web app](https://img.shields.io/badge/web%20app-gmg.fungalgenomes.org-blue)](https://gmg.fungalgenomes.org/)
[![Container](https://img.shields.io/badge/container-ghcr.io%2Fstajichlab%2Fgmg-blue)](https://github.com/stajichlab/GmG/pkgs/container/gmg)
[![License: MIT](https://img.shields.io/github/license/stajichlab/GmG)](LICENSE)

Given a fungal genome (or reads), guess where it belongs: phylum, class, order, family, genus, species, plus an
ANI-cluster, each with a confidence. Built for MATPredict routing, meant as a general tool.

Status: working CLI and core (`gmg`), a WASM build, and a 23,669-genome BFD reference built on the HPCC.
Results: [reports/overnight_2026-10-27.md](reports/overnight_2026-10-27.md). History: [LOG.md](LOG.md). Design: [plan.md](plan.md).

## Install

Pick one. All of them run the same code.

**Python** (Linux x86_64 / aarch64 and macOS Apple Silicon, Python 3.9 or newer; the import name is `gmg`):

```
pip install guess-my-genome
python -c "import gmg; print(gmg.__version__)"
```

`pip install gmg` is a different, unrelated project. Other platforms (Intel Mac, Windows): build from this repository with maturin ([crates/gmg-py](crates/gmg-py/README.md)) or use a binary or the container.

**Command-line binary** (Linux x86_64 and macOS Apple Silicon), from the [latest release](https://github.com/stajichlab/GmG/releases/latest):

```
V=v0.3.1; P=osx-arm64                  # or P=linux-64
curl -fsSLO https://github.com/stajichlab/GmG/releases/download/$V/gmg-$V-$P.tar.gz
curl -fsSLO https://github.com/stajichlab/GmG/releases/download/$V/gmg-$V-$P.tar.gz.sha256
shasum -a 256 -c gmg-$V-$P.tar.gz.sha256     # on Linux: sha256sum -c
tar xzf gmg-$V-$P.tar.gz && ./gmg-$V-$P/gmg --version
```

**Container** (no install; `:X.Y.Z` pins a release, `:latest` is the newest release, `:edge` is the tip of `main`):

```
docker run --rm -v "$PWD:/data" ghcr.io/stajichlab/gmg:latest run /data/genomes --index /data/bfd.gmgi -o /data/results
```

**In the browser:** <https://gmg.fungalgenomes.org/> runs the same code as WebAssembly on your computer; nothing is uploaded ([docs/web.md](docs/web.md)).

**A reference index** is needed for `run` and `guess`. Indexes are released separately from the code
([index-distribution.md](docs/index-distribution.md)): `bfd-reps-*-s20000.gmgi` (about 67 MB, fast, fewer species calls) or the full `bfd-*-s10000.gmgi` (about 213 MB):

```
gh release download index-2026.10.5 -R stajichlab/GmG --pattern 'bfd-reps-*.gmgi'     # or download from the release page
```

From source (no Homebrew or rustup needed): see [Develop](#develop-no-homebrew-or-rustup-needed).

## Quick start

```
gmg run my_genomes/ --index bfd-reps-2026.10.5-s20000.gmgi -o results/     # a folder (or files) of FASTA / FASTQ
column -t -s$'\t' results/results.tsv | less -S                              # one row per input; hits in results/hits.tsv
```

`results/report.html` is the report (open it in a browser); `results.tsv` has one row per input with the calls from phylum to species, confidence, QC flags and status; `run.json` records how the
run was made. Columns are documented in [docs/results-schema.md](docs/results-schema.md). Other commands: `gmg sketch / build / guess / info / convert / cluster / bench / completions / man`.
In the browser: [docs/web.md](docs/web.md). Example genomes to try: [examples/](examples/README.md). From Python: [crates/gmg-py](crates/gmg-py/README.md). Step by step: [docs/tutorial.md](docs/tutorial.md); what the results do and do not mean: [docs/limitations.md](docs/limitations.md); every command has examples in `gmg <command> --help`.

## A set of read files: `popsketch`

`popsketch run reads/ -o out/ --labels labels.tsv` sketches every sample (R1/R2 merged), audits the set against itself (samples that are another organism, low coverage, a second genome mixed in, half the content missing, a label the data contradict), and for the samples that pass writes a PCA, clusters, a first-pass tree and an HTML report. It shows coarse structure, not a phylogeny, and it states where the data stop being informative. See [docs/popsketch.md](docs/popsketch.md) and the [tutorial](docs/popsketch-tutorial.md) (`popsketch demo` writes a small example); the design and its evidence are in [docs/popstruct-plan.md](docs/popstruct-plan.md).

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
| `crates/gmg-input` | finding FASTA/FASTQ files, naming samples, merging R1/R2 and lanes, sample sheets (shared by `gmg run` and `popsketch`) |
| `crates/popsketch` | `popsketch run`: audit a set of read files and show coarse population structure ([docs/popsketch.md](docs/popsketch.md)) |
| `crates/gmg-report` | report engine: standalone HTML report with inline SVG charts (no I/O; also runs in WASM) |
| `crates/gmg-wasm` | JavaScript API (`Index`: load, sketch, guess; `render_report`) |
| `scripts/` | analysis and audit (Python), `scripts/hpcc/` SLURM jobs |
| `reports/`, `data/` | results, calibrated thresholds, cluster tables, label audit |

## Using GmG from MATPredict

[MATPredict](https://github.com/stajichlab/MATPredict) asks GmG which phylum a genome is in so it only searches that phylum's curated mating-type loci:
`matpredict detect --genome G.fna --out-dir out/ --phylum auto` (and `auto-weak`, which also accepts a weak call; the browser tool has a "Guess phylum" button). Every
record has a `routing_hint` (`usable`, `phylum`, the deepest called taxid; [docs/results-schema.md](docs/results-schema.md)) that says whether the phylum call is safe to route on:
it is not usable when there is no call, the input looks non-fungal, or the call is weak (best-hit ANI below 0.85). MATPredict needs the `gmg` binary or the Python bindings.

## HPCC

A checkout lives in `~/projects/GmG` (`gh repo clone stajichlab/GmG`); work is committed here, pushed, then `git pull` there.
Build with `pixi run build` in that checkout (login node, short compile) and install with `scripts/hpcc/deploy_gmg.sh`
(versioned binary + symlink swap; never overwrite a running binary). Data and results are under
`/bigdata/stajichlab/shared/projects/GmG_runs/`.

## Releases and container

Tagging `vX.Y.Z` runs two workflows: `release.yml` (`gmg` binaries for linux-64 and osx-arm64, the wasm package and an OCI image,
`ghcr.io/stajichlab/gmg:X.Y.Z` and `:latest`) and `pypi.yml` (wheels, published to [PyPI](https://pypi.org/project/guess-my-genome/) with trusted publishing after
approval of the `pypi` environment). Every merge to `main` also publishes `:edge`. Maintainers: [docs/releasing.md](docs/releasing.md). Image:

```
docker run --rm -v "$PWD:/data" ghcr.io/stajichlab/gmg:edge guess --index /data/bfd.gmgi /data/genome.fa
apptainer pull gmg.sif docker://ghcr.io/stajichlab/gmg:edge      # on the HPCC (module load apptainer, on a compute node)
```
