#!/usr/bin/env python3
"""Truth table and file list for the Coccidioides read validation (from MATPredict's manifest of the 2025_All_Cocci read sets).

usage: make_cocci_truth.py MANIFEST.tsv OUT_PREFIX
writes OUT_PREFIX.truth.tsv (id = R1 file name without .fastq.gz; lineage; strain; relabelled) and OUT_PREFIX.r1.list (R1 paths on the HPCC).
Species are the CORRECTED species of the manifest (two strains were relabelled); `relabelled` = yes marks them.
"""
import csv, sys

base = "/bigdata/stajichlab/shared/projects/Coccidioides/PopGenomics/2025_All_Cocci/Assembly/input/"
rows = [r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t") if r["read1"]]
with open(sys.argv[2] + ".truth.tsv", "w") as t, open(sys.argv[2] + ".r1.list", "w") as l:
    t.write("id\tphylum\tclass\torder\tfamily\tgenus\tspecies\tstrain\trelabelled\tbytes\n")
    for r in rows:
        rid = r["read1"].replace(".fastq.gz", "")
        t.write(f"{rid}\tAscomycota\tEurotiomycetes\tOnygenales\tOnygenaceae\tCoccidioides\t{r['species']}\t{r['strain']}\t{r['relabelled']}\t{r['bytes']}\n")
        l.write(base + r["read1"] + "\n")
print(len(rows), "read sets;", sum(1 for r in rows if r["species"] == "Coccidioides immitis"), "immitis,", sum(1 for r in rows if r["species"] == "Coccidioides posadasii"), "posadasii")
