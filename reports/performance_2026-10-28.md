# Performance review (2026-10-28)

Question: is algorithmic work needed for speed? Measured with the release build (macOS arm64, 8 cores) on the published
full index (`bfd-2026.10-s10000.gmgi`, 23,669 genomes, 199 MB) unless noted.

| what | measurement | verdict |
|---|---|---|
| Sketching assemblies, native | 115.6 MB of genomes (5 fungi) in 1.0 s of CPU after index load: ~130 MB/s per thread | linear, fast; the BFD collection (~23.7k genomes, ~0.7 TB) is ~1.5 CPU-hours |
| Sketching reads, native, gz | 2M x 150 bp reads (113 MB gz, 300 Mb) in 3.9 s on 1 thread: ~85 Mb/s (uncompressed: 115 Mb/s) | fast |
| Search + calling | 13 ms for a 1,100-hash query against the full index (WASM); invisible natively | inverted index touches only the query's hashes |
| Index load, native | 0.9 s, 776 MB resident for the 199 MB file; representatives index (61 MB): 0.3 s | acceptable for batch runs; matters for many short invocations |
| Index load, WASM (node) | 0.3 s / 0.9 s; **peak resident 450 MB (reps) and 2.1 GB (full)** | the real limit for the browser |
| Sketching, WASM | 56-72 MB/s (about half of native) | fine: a 30 Mb genome in under a second |
| 1,257 assemblies on the HPCC | ~95 s for 194 stress inputs on 4 threads; the 51 independent-set batches ran in ~13 s each | fine |
| Reads on the HPCC | 558 sets x 2M reads in 1,896 s on 8 threads = ~27 CPU-s per set = ~11 Mb/s per thread, **8x slower than the same work locally** | I/O (network file system, 760 MB gz files read from the start) or CPU contention, not the algorithm; not investigated |

## Conclusion

No algorithmic optimisation is needed for classification: sketching is linear in the input and the search cost is
proportional to the query's hashes, both measured far below anything users will notice. Remaining work is engineering:

1. **Memory and cold start of the index (issue #36).** The 199 MB file is expanded to ~560 MB of arrays and the loader
   first reads the whole file into a second buffer (~2 GB peak in the browser). Options, cheapest first: decode
   straight into the arrays without the intermediate buffer; keep postings delta-coded and decode a list only when a
   query hash hits it (a hash directory plus byte offsets lets a native build `mmap` the file and load in milliseconds, and
   a browser build touch only the pages it needs); ship the representatives index (450 MB peak) for the browser.
2. **Reads on shared storage:** check whether the HPCC slowdown is I/O (copy a read file to local scratch and compare)
   before touching code.
3. Possible later: faster gzip (libdeflate / zlib-ng backend, 2-3x on reads) and splitting one huge FASTA across threads;
   neither is needed at current sizes.

## Where algorithmic work will matter (not started)

* PopSketch strain trees / PCA (#15): all-vs-all containment is O(N^2) pairs; N = 15 is trivial, N = 10^4 needs a sparse
  inverted-index pass instead of pairwise merges, a neighbour-joining variant below O(N^3) and a randomized SVD for PCA.
* Mixed-sample decomposition (#31): greedy set cover of the query's hashes over reference genomes (sourmash gather).
  With the inverted index each greedy step is incremental, so the cost is manageable, but the design (abundance
  weighting, error-prone hashes, stopping rule) is an algorithm-design task.

These are the places where a deeper model for design review is worth it; classification speed is not.
