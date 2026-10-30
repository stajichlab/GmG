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
gmg run examples_data/ --index bfd-reps-2026.10.2-s20000.gmgi -o example_results/
```

Expected: the five fungi are called to genus or species with their own phylum (Ascomycota x2, Basidiomycota, Wallemiomycetes
under Basidiomycota, Microsporidia); *E. coli* is `no_match`. If your index has a different taxonomy version the names may differ
slightly (see docs/taxonomy.md).
