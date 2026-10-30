#!/usr/bin/env python3
"""Give the clusters whose centroid was removed from the reference a new centroid (#47, #48).

usage: reseed_clusters.py META_IN.tsv META_OUT.tsv EXCLUDED.txt [ORDER.txt]
META_IN is a meta_full.tsv that still carries the `cluster` column (centroid id) of every genome, META_OUT the same table without the excluded
genomes and with `cluster` rewritten where the centroid is gone: the new centroid is the first remaining member in ORDER.txt (the BFD priority
order the clustering used, index/reps.order.txt; default: the order of the table). Members keep their cluster; only the label changes. Prints
what changed. Run it in `15_release_index.sbatch` in place of the plain exclusion filter, otherwise `cluster` points at a genome that is not in the index.
"""
import csv, sys

meta_in, meta_out, excl_f = sys.argv[1:4]
order_f = sys.argv[4] if len(sys.argv) > 4 else None
ex = {l.strip() for l in open(excl_f) if l.strip()}
rows = list(csv.reader(open(meta_in), delimiter="\t"))
head, body = rows[0], rows[1:]
ic, ii = head.index("cluster"), head.index("id")
keep = [r for r in body if r[ii] not in ex]
rank = {}
if order_f:
    for n, l in enumerate(open(order_f)):
        rank.setdefault(l.split("\t")[0].strip(), n)
members = {}
for r in keep:
    members.setdefault(r[ic], []).append(r[ii])
new = {}
for c, ms in members.items():
    if c in ex:
        new[c] = min(ms, key=lambda i: rank.get(i, 10**9))
for r in keep:
    if r[ic] in new:
        r[ic] = new[r[ic]]
with open(meta_out, "w") as f:
    f.write("\t".join(head) + "\n")
    for r in keep:
        f.write("\t".join(r) + "\n")
print(f"excluded {len(body) - len(keep)} genomes; {len(new)} clusters re-seeded ({sum(len(members[c]) for c in new)} members): "
      + ", ".join(f"{c[:20]}->{n[:20]}" for c, n in list(new.items())[:5]), file=sys.stderr)
