#!/usr/bin/env python3
"""D1h: do collapsed clusters mix labels? (A label test on representatives gives a cluster its representative's label.)
usage: label_purity.py SLUG collapse.tsv labels.tsv >> label_purity.tsv
Columns: slug column clusters multi_member mixed_clusters mixed_share_of_multi samples_in_mixed
 multi_member: clusters with at least 2 samples; mixed: of those, clusters whose members carry at least 2 different non-empty labels
 samples_in_mixed: samples whose label differs from their representative's label."""
import csv
import sys
from collections import defaultdict


def main():
    slug, col_path, lab_path = sys.argv[1:4]
    lab = {}
    with open(lab_path) as f:
        r = csv.DictReader(f, delimiter="\t")
        cols = [c for c in r.fieldnames if c != "id"]
        for row in r:
            lab[row["id"]] = row
    members = defaultdict(list)
    with open(col_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            members[row["representative"]].append(row["id"])
    for c in cols:
        multi = [m for m in members.values() if len(m) >= 2]
        mixed = [m for m in multi if len({lab[i][c] for i in m if lab.get(i, {}).get(c)}) >= 2]
        differ = sum(1 for rep, m in members.items() for i in m
                     if lab.get(i, {}).get(c) and lab.get(rep, {}).get(c) and lab[i][c] != lab[rep][c])
        print(f"{slug}\t{c}\t{len(members)}\t{len(multi)}\t{len(mixed)}\t{len(mixed)/max(len(multi),1):.2f}\t{differ}", flush=True)


if __name__ == "__main__":
    main()
