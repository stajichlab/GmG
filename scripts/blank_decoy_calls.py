#!/usr/bin/env python3
"""Treat a phylum call that only says "nearest to a non-fungal decoy, but distantly" as no call (scoring helper, #39).

usage: blank_decoy_calls.py RESULTS.tsv INDEX_META.tsv DECOYS.tsv OUT.tsv
Rows with the QC flag `distant_non_fungal_hit` whose phylum call is one of the decoy phyla (the phyla of the genomes listed in
DECOYS.tsv, as the index holds them) get every call_* column blanked and `phylum_call_decoy` set: the flag's own documentation says
"not a verdict; treat the call as unknown" (docs/qc.md). Everything else is unchanged. Prints how many rows and what truth they had.
"""
import csv, sys

res_p, meta_p, dec_p, out = sys.argv[1:5]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
dec_ids = {r["id"] for r in csv.DictReader(open(dec_p), delimiter="\t")}
phyla = {r["phylum"] for r in csv.DictReader(open(meta_p), delimiter="\t") if r["id"] in dec_ids}
rd = csv.DictReader(open(res_p), delimiter="\t")
rows = list(rd)
n = 0
for r in rows:
    r["phylum_call_decoy"] = ""
    if "distant_non_fungal_hit" in r["qc_flags"] and r["call_phylum"] in phyla:
        r["phylum_call_decoy"] = r["call_phylum"]
        for k in RANKS:
            r["call_" + k] = ""
        n += 1
with open(out, "w") as o:
    w = csv.DictWriter(o, fieldnames=rd.fieldnames + ["phylum_call_decoy"], delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(f"{res_p}: {n} of {len(rows)} calls blanked (distant decoy hit; {len(phyla)} decoy phyla)")
