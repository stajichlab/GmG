#!/usr/bin/env python3
"""Sample sheet for `gmg run` on the truth-set queries (tool comparison, #46): make_sheet.py S2_TRUTH.tsv QUERY_DIR MODE (loo|los|log).
Columns id, path, expected (the truth species or genus, used to hide the species/genus when the genome is not in the index), ref (the BFD
genome this query is a copy of)."""
import csv, sys

truth, qdir, mode = sys.argv[1:4]
print("id\tpath\texpected\tref")
for r in csv.DictReader(open(truth), delimiter="\t"):
    exp = r["species"] if mode == "los" else r["genus"] if mode == "log" else r["species"]
    print(f"{r['id']}\t{qdir}/{r['id']}.fa.gz\t{exp}\t{r['bfd_id']}")
