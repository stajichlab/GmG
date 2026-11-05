#!/usr/bin/env python3
"""samples.csv -> manifest.tsv (id, organism, notes, R1, R2, found). Run on the HPCC."""
import csv, os, re, sys
base = sys.argv[1]  # .../2025_All_Cocci/Assembly
out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
out.writerow(["id", "organism", "notes", "r1", "r2", "found"])
for r in csv.DictReader(open(os.path.join(base, "samples.csv"))):
    fb = r["FileBase"] or "_R[12]_001.fastq.gz"
    acc = r["RunAcc"]
    m = re.search(r"SRR\d+", acc)  # e.g. NACVFR_CocciSRR21204713 -> SRR21204713
    cands = [acc] + ([m.group(0)] if m else [])
    for a in cands:
        r1 = os.path.join(base, "input", a + fb.replace("[12]", "1"))
        if os.path.exists(r1):
            break
    r2 = os.path.join(base, "input", a + fb.replace("[12]", "2"))
    out.writerow([r["Strain"], r["Organism"], r["Notes"], r1, r2, int(os.path.exists(r1))])
