#!/usr/bin/env python3
"""Pick the 12 smallest runs of 7-20 Gb for each species into depth/tasks.tsv (slug, run, base_count). Run from GmG_popgen."""
import csv, os, sys
os.makedirs("depth", exist_ok=True)
with open("depth/tasks.tsv", "w") as o:
    for slug in sys.argv[1:]:
        rows = [r for r in csv.DictReader(open(f"popsra/{slug}/runs.tsv"), delimiter="\t") if 7e9 <= int(r["base_count"]) <= 20e9]
        rows.sort(key=lambda r: int(r["base_count"]))
        for r in rows[:12]:
            o.write(f'{slug}\t{r["run_accession"]}\t{r["base_count"]}\n')
            print(slug, r["run_accession"], round(int(r["base_count"]) / 1e9, 1), r["study_accession"])
