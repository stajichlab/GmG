#!/usr/bin/env python3
"""Join BFD samples.csv to the genome files and write the inputs for `gmg sketch` / `gmg build`.

usage: make_meta.py SAMPLES_CSV GENOME_DIR OUT_DIR

writes OUT_DIR/meta.tsv (id taxid phylum class order family genus species),
       OUT_DIR/genomes.list (absolute .fa.gz paths that have a samples.csv row),
       OUT_DIR/files_without_row.txt, OUT_DIR/rows_without_file.txt
"""
import csv, os, sys

samples, gdir, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
files = {f[:-len(".fa.gz")]: f for f in os.listdir(gdir) if f.endswith(".fa.gz")}
rows = {}
with open(samples, newline="") as fh:
    for r in csv.DictReader(fh):
        rows[r["ASMID"]] = r

both = sorted(set(files) & set(rows))
with open(os.path.join(out, "meta.tsv"), "w") as m, open(os.path.join(out, "genomes.list"), "w") as g:
    m.write("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
    for i in both:
        r = rows[i]
        tid = r["NCBI_TAXONID"] or "0"
        f = [r[c].strip().replace("\t", " ") for c in ("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")]
        m.write("\t".join([i, tid] + f) + "\n")
        g.write(os.path.join(gdir, files[i]) + "\n")
with open(os.path.join(out, "files_without_row.txt"), "w") as f:
    f.write("\n".join(sorted(set(files) - set(rows))) + "\n")
with open(os.path.join(out, "rows_without_file.txt"), "w") as f:
    f.write("\n".join(sorted(set(rows) - set(files))) + "\n")
print(f"files={len(files)} rows={len(rows)} matched={len(both)} "
      f"files_without_row={len(set(files)-set(rows))} rows_without_file={len(set(rows)-set(files))}")
