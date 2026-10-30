# Performance, refreshed for index 2026.10.2 (2026-10-30, #49)

Same questions as `performance_2026-10-28.md`, measured with `scripts/perf_native.sh` on the published indexes of release `index-2026.10.2`
(full: 23,750 genomes, 213 MB; representatives: 7,529 genomes, 67 MB; both with the 81 decoy genomes) and the 17-genome golden panel
(274 MB of sequence in gzipped FASTA, 90 MB on disk). One run per cell; wall time and peak resident memory (`/usr/bin/time`); each `gmg run` includes
reading and decompressing the files, sketching, searching, calling and writing the TSVs (no report).

| what | macOS arm64 (Apple M3, 8 cores) | Linux x86-64 (HPCC node r41, 8 cpus) |
|---|---|---|
| `gmg info`, full index (load) | 0.20 s, 226 MB | 0.27 s, 228 MB |
| `gmg info`, representatives index | 0.06 s, 75 MB | 0.14 s, 75 MB |
| classify the panel, full index, 1 thread | 3.5 s, 288 MB | 6.2 s, 250 MB |
| classify the panel, representatives index, 1 thread | 3.3 s, 136 MB | 6.0 s, 101 MB |
| classify the panel, full index, 8 threads | 1.0 s, 381 MB | 1.9 s, 339 MB |
| classify the panel, representatives index, 8 threads | 0.8 s, 230 MB | 1.7 s, 190 MB |
| sketch only (`gmg sketch`, k = 21, scaled 10000), 1 thread | 3.3 s, 61 MB | 5.9 s, 32 MB |

* **Throughput:** about 80 MB of gzipped-FASTA sequence per second per thread on the Mac and 46 MB/s on the cluster node, including decompression;
  sketching is almost the whole cost (sketch only 3.3 s of 3.5 s), the search and the calling add 0.2 s for 17 genomes. 8 threads give 3.5x (Mac) and 3.2x (Linux).
  At the 8-thread rate (141-277 MB/s) the 1,257 post-freeze assemblies (roughly 40 GB of sequence, estimated from ~35 Mb per genome) would take 3-5 minutes of compute on one
  8-core node, plus the download; this was not timed as one run.
* **Memory:** the full index costs ~230 MB at load (the packed v3 format stays near the file size, as in the 2026-10-29 update), a classification run
  250-380 MB, the representatives index 75-230 MB. A laptop and the browser are fine with the representatives index; the full index fits any node.
* **Against 2026-10-28/29** (full index, same format): load 0.09 s / 220 MB then, 0.20 s / 226 MB now on a busier machine; the index grew 7% (213 vs 199 MB: decoys, taxids) and the memory
  with it by about 3%. No regression in per-genome cost was found; the cells are single runs and the workstation was running other builds, so differences below ~30% are noise.
* **Not re-measured here:** WebAssembly peak memory and speed (the 2026-10-29 numbers (463 MB for the full index, 178 MB for the representatives index) were for the
  index without decoys and are not repeated for this one), and reads (a set of 2M reads, 73 MB with the representatives index; the HPCC reads slowdown is storage, see the first report).
  The WASM build is checked for identical results on the golden panel (`reports/validation/reproducibility.md`).

Reproduce: `scripts/perf_native.sh GMG FULL.gmgi REPS.gmgi PANEL_DIR` (the panel is cached by `python3 scripts/golden_check.py --fetch-only --cache DIR`);
on the cluster `scripts/hpcc/29_perf_linux.sbatch`.
