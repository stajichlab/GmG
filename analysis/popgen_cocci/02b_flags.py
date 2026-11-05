#!/usr/bin/env python3
"""02b_flags.py QC.tsv > flagged.txt : the 9 reported strains + samples whose size at decent coverage is off (>32 Mb)."""
import csv, sys
rep = "B10918 B10919 B11035 B11036 B11057 B11198 B11871 B12399 B12526".split()
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    if "Coccidioides" in r["organism"] and (r["id"] in rep or (float(r["est_genome_mb"]) > 32 and float(r["median_kmer_cov"]) >= 5)):
        print(r["id"])
