# GmG: Guess My Genome

Given a fungal genome (later: reads), guess where it belongs: phylum, order, genus, species,
with a confidence for each rank. Built for MATPredict routing but meant as a general tool.

Status: early. `crates/gmg-core` is a sourmash-compatible FracMinHash sketching and comparison
core (no I/O, wasm-friendly). Everything else is in [plan.md](plan.md).

```
cargo test --release            # includes a hash-for-hash check against the sourmash crate
                                # (needs data/genomes/, see data/accessions.txt; skipped if absent)
```

On macOS with Homebrew Rust, wasm builds need the rustup toolchain first on PATH:
`PATH=$HOME/.rustup/toolchains/stable-aarch64-apple-darwin/bin:$PATH cargo build --release --target wasm32-unknown-unknown -p gmg-core`
