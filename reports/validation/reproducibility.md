# Reproducibility on real genomes (#49)

Does the same genome give the same hashes and the same calls on every platform and against an independent FracMinHash implementation? The test
set is the **golden panel**: 17 real assemblies (`data/golden/panel.tsv`: five fungal lineages, a bacterium, a decoy, a low-confidence and a divergent
input; each file pinned by SHA-256) with the calls, flags and ANI they must give (`data/golden/expected.tsv`) for the pinned published index
(`bfd-reps-2026.10.2-s20000.gmgi`, `data/golden/INDEX.txt`). Check: `scripts/golden_check.py` (also a CI job).

| leg | what was compared | result |
|---|---|---|
| macOS arm64 (Apple M3), release build | `gmg run` on the panel against `expected.tsv` | 17 / 17 identical (the file was made here) |
| Linux x86-64, GitHub Actions (ubuntu-latest) | the same, binary built from the commit under test | **17 / 17 identical** (`golden.yml`, every change under `crates/`, weekly) |
| macOS arm64, GitHub Actions (macos-14) | the same | **17 / 17 identical** |
| Linux x86-64, UCR HPCC, the deployed binary | the same | **17 / 17 identical** (`scripts/hpcc/28_parity_linux.sbatch`) |
| WebAssembly (Node), current build | records of `Index.classify` against the records of `gmg run` for all 17 genomes, every field, numbers to 1e-12 relative | **17 / 17 equal** (`scripts/wasm_golden.mjs`; the CI parity test covers the synthetic inputs) |
| Python sourmash 4.9.3 (HPCC module), k = 21, scaled = 10000, seed 42 | the sorted hash set of `gmg sketch` against `sourmash sketch dna` for each genome | **17 / 17 identical hash sets** (424 to 3,394 hashes per genome; `scripts/sourmash_parity.py`) |
| Container (ghcr.io/stajichlab/gmg) | not run here: no container runtime on the HPCC and none running on the workstation | CI builds the image and smoke-tests it; the binary inside is the Linux build above. Not independently verified on the panel |

What each leg shows: the sketch (hash set) is bit-identical to an independent implementation (sourmash); the classification gives identical status, calls at
every rank, QC flags, top hit and ANI (to 0.0005 in the golden check; to 1e-12 relative in the wasm comparison) on two operating systems, two CPU
architectures and wasm. The tolerances exist because libm may differ in the last bit between wasm and native; no difference larger than the tolerance was found.

## What the golden check caught

The first CI run **failed on both platforms**: one panel row (*E. coli*, a `no_match`) differed in flag, top hit and ANI. Cause: `distant_non_fungal_hit`
had been changed on `main` (#93, only flag when a phylum call was made) after the expectations were generated; the check did its job. The
expectation was regenerated from the new code, after confirming that the one row changed, and the three fields that changed are
exactly those that #93 should change. Changing an expectation is deliberate: `--update`, with the reason in the commit message.

## Notes

* sourmash 4.9.3 is the version on the HPCC module system; the roadmap (#16) names 4.9.4. The sketching code of 4.9.3 and 4.9.4 is the same (a Python
  3.10+ environment is needed for 4.9.4; the workstation has 3.9). Re-check with 4.9.4 when the CI image has it.
* The panel pins the **code** (calls must not change for a given index). When a new index is published the expectations are re-pinned with `--update`
  and the diff is read (see `docs/recalibration.md`).
* Index loading is not part of the comparison (all legs read the same published file and verify its checksum).
