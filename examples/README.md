# Example genomes

Six small public assemblies to try GmG with: five fungi from different lineages and one bacterium (a negative control).
They are listed in `examples.tsv` (accession, name, what to expect, NCBI FTP URL); nothing is stored in this repository.
Five of the six are recent assemblies that are **not** in the reference, so the calls are genuine classifications (not
`identical_to_reference`). Total download: about 17 MB.

* **Web app:** under "Genomes or reads", the *Try an example* buttons download the file from NCBI into your browser and add it
  to the queue (NCBI's FTP allows cross-origin downloads; nothing is sent anywhere).
* **Command line:**

```
scripts/fetch_examples.sh examples_data          # downloads all six (or: scripts/fetch_examples.sh DIR scer ecoli)
gmg run examples_data/ --index bfd-reps-2026.10.3-s20000.gmgi -o example_results/
```

Measured with `gmg` from `main` and the 2026.10.3 indexes (a later index or taxonomy may differ slightly, see docs/taxonomy.md):

| example | representatives index (67 MB, the web app) | full index (213 MB) |
|---|---|---|
| *Saccharomyces cerevisiae* | genus *Saccharomyces* (species withheld: the nearest representatives disagree on the species) | *Saccharomyces cerevisiae*, ANI 0.9985 |
| *Schizosaccharomyces pombe* | *S. pombe*, ANI 0.9943 | *S. pombe*, ANI 1.0000 |
| *Malassezia restricta* | genus *Malassezia* | *M. restricta*, ANI 0.9997 |
| *Wallemia ichthyophaga* | *W. ichthyophaga*, ANI 0.9963 | same, 0.9996 |
| *Vavraia culicis* | *V. culicis*, ANI 0.9979 | same, 0.9980 |
| *Escherichia coli* K-12 | `no_match` (ANI 0.80, 2 shared hashes: not a fungus) | `no_match` |

The yeast and *Malassezia* rows show what the representatives index gives up: only the best genome per species is kept, so a species with several
close relatives in the reference can fall back to the genus. Use the full index when the species matters.
