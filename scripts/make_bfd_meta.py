#!/usr/bin/env python3
"""Write data/external/bfd_meta.tsv (BFD label of every reference genome) from BFD's tables/species.parquet.

usage: make_bfd_meta.py SPECIES_PARQUET OUT_TSV

Columns match scripts/make_meta.py (id taxid phylum class order family genus species), sorted by id, so the
file can be diffed against the previous one. SUBPHYLUM and SUBCLASS are in the parquet but are not written:
select_external.py, select_post_freeze.py and make_upstream_report.py read this file by position.
Run on the HPCC (needs pyarrow or duckdb).
"""
import csv, sys

src, out = sys.argv[1:3]
try:
    import pyarrow.parquet as pq
    rows = pq.read_table(src).to_pylist()
except ImportError:
    import duckdb
    cur = duckdb.connect().execute("select * from read_parquet(?)", [src])
    cols = [d[0] for d in cur.description]
    rows = [dict(zip(cols, r)) for r in cur.fetchall()]

def clean(v):
    return ("" if v is None else str(v)).strip().replace("\t", " ")

seen = {}
for r in rows:
    seen[clean(r["ASMID"])] = r
with open(out, "w") as fh:
    fh.write("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
    for i in sorted(seen):
        r = seen[i]
        tid = clean(r["NCBI_TAXONID"])
        tid = tid.split(".")[0] if tid else "0"
        f = [clean(r[c]) for c in ("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")]
        fh.write("\t".join([i, tid or "0"] + f) + "\n")
print(f"rows={len(rows)} written={len(seen)}")
