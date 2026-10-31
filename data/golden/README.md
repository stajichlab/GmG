# Golden panel (#49)

17 real genomes (`panel.tsv`: accession, NCBI FTP URL, SHA-256 of the file; an assembly version never changes) and the calls `gmg` must give for them
(`expected.tsv`) with the pinned index (`INDEX.txt`): status, best rank, the call at every rank, QC flags, top hit and ANI (to 0.0005).

```
python3 scripts/golden_check.py --gmg target/release/gmg --index bfd-reps-2026.10.5-s20000.gmgi     # exit 1 and a diff on any change
python3 scripts/golden_check.py ... --update       # rewrite expected.tsv: only for an intended change, with the reason in the commit message
```

The panel mixes lineages (Ascomycota yeasts and filamentous fungi, Basidiomycota, Mucoromycota, Microsporidia), a bacterium (no_match), a decoy
(`not_fungal`), a low-confidence input and a divergent assembly with no close relative. `.github/workflows/golden.yml` runs the check on Linux and macOS for
every change under `crates/` and weekly; `reports/validation/reproducibility.md` has the cross-platform, WASM and sourmash results.
