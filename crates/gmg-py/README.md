# gmg: Python bindings (#37)

The same Rust code as `gmg run` and the web app, for notebooks and scripts.

```python
import gmg
ix = gmg.Index("bfd-reps-2026.10.3-s20000.gmgi")     # any index from the releases
rec = ix.classify("genome.fa.gz")                     # a dict: the record `gmg run` writes to results.jsonl
rec["status"], rec["best_rank"]                       # ok | low_confidence | no_match | error
rec["guess"]["calls"]["genus"]["taxon"]               # per-rank calls, ANI, top hits, QC flags in rec["qc"]

records = [ix.classify(p) for p in paths]             # failures are records (status "error"), not exceptions
open("results.tsv", "w").write(gmg.results_tsv_text(records))      # identical to the CLI's results.tsv

h = ix.sketch(open("x.fa", "rb").read())              # sorted unique hashes at the index's k / scaled
ix.guess_hashes(h)                                    # the JSON of `gmg guess`, as a dict
ix.use_thresholds("novel")                            # known (default) | permissive | novel
```

`Index.classify(path, kind="auto", min_abund=0, max_records=0, top=5, id=None)` reads plain or gzipped FASTA/FASTQ; `classify_bytes(data, name=...)`
takes sequence text in memory. Not in the bindings: sample sheets, read pairing and `--exclude-*` (use the CLI for those).

## Build and test

```
pip install maturin pytest
cd crates/gmg-py && maturin develop          # or: maturin build --release (an abi3 wheel for Python >= 3.9)
GMG_BIN=../../target/debug/gmg pytest tests  # the tests build a small index with the gmg binary and compare with `gmg run`
```

The crate is outside the main Cargo workspace (like `fuzz/`), so `cargo test` at the root does not need Python; `.github/workflows/python.yml` builds and tests it.
