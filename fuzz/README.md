# Fuzz targets (#34)

`index` (any bytes as an index file: streamed and packed loaders, then a search), `sketch_file` (any bytes as a sketch file) and
`fastx` (any bytes as FASTA/FASTQ, the entry points the web app calls). The pass criterion is "an error or a result, never a
panic, abort or huge allocation". The deterministic version of these checks (mutated and truncated files, 25 hostile texts)
runs in every CI pass (`crates/gmg-core/tests/robust.rs`); this adds coverage-guided search.

```
rustup toolchain install nightly && cargo install cargo-fuzz
cargo +nightly fuzz run index -- -max_total_time=60      # also: sketch_file, fastx
```

CI (`.github/workflows/fuzz.yml`) runs each target for 60 s on pull requests that touch `crates/gmg-core` and nightly for 10 min;
a crash uploads the input as an artifact. Reproduce with `cargo +nightly fuzz run TARGET path/to/crash-file`.
