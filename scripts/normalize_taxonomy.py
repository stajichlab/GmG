#!/usr/bin/env python3
"""Rebuild the lineage of every reference genome from its NCBI taxid (issue #30).

usage: normalize_taxonomy.py META.tsv TAXDUMP_DIR OUT.tsv [--report REPORT.md] [--clusters CLUSTERS.tsv] [--taxdump-label TEXT]

META.tsv      id taxid phylum class order family genus species   (BFD labels, from make_meta.py)
TAXDUMP_DIR   directory with nodes.dmp names.dmp merged.dmp (NCBI new_taxdump)
OUT.tsv       id taxid phylum class order family genus species cluster phylum_taxid ... species_taxid label taxid_status
              ranks and names come from NCBI; `label` is the original BFD species string; cluster comes from --clusters
Taxids that were merged are followed; a taxid that is deleted or unknown keeps its BFD labels (taxid_status says so).
"""
import argparse, csv, collections, sys, os

RANKS = ["phylum", "class", "order", "family", "genus", "species"]
ap = argparse.ArgumentParser()
ap.add_argument("meta"); ap.add_argument("taxdump"); ap.add_argument("out")
ap.add_argument("--report"); ap.add_argument("--clusters"); ap.add_argument("--taxdump-label", default="")
a = ap.parse_args()

parent, rank, name, merged, deleted = {}, {}, {}, {}, set()
with open(os.path.join(a.taxdump, "nodes.dmp")) as f:
    for l in f:
        p = l.split("\t|\t")
        parent[int(p[0])] = int(p[1]); rank[int(p[0])] = p[2]
with open(os.path.join(a.taxdump, "names.dmp")) as f:
    for l in f:
        p = l.rstrip("\n").split("\t|\t")
        if p[3].startswith("scientific name"):
            name[int(p[0])] = p[1]
with open(os.path.join(a.taxdump, "merged.dmp")) as f:
    for l in f:
        p = l.split("\t|\t"); merged[int(p[0])] = int(p[1].strip("\t|\n "))
dpath = os.path.join(a.taxdump, "delnodes.dmp")
if os.path.exists(dpath):
    with open(dpath) as f:
        deleted = {int(l.split("\t|")[0]) for l in f}

clusters = {}
if a.clusters:
    for r in csv.DictReader(open(a.clusters), delimiter="\t"):
        clusters[r["id"]] = r["cluster"]


def lineage(t):
    """({rank: (taxid, name)}, status)"""
    seen = 0
    while t in merged and seen < 10:
        t = merged[t]; seen += 1
    if t not in parent:
        return {}, ("deleted" if t in deleted else "unknown")
    out = {}
    cur = t
    for _ in range(100):
        r = rank.get(cur)
        if r in RANKS and r not in out:
            out[r] = (cur, name.get(cur, ""))
        if parent[cur] == cur:
            break
        cur = parent[cur]
    return out, ("merged" if seen else "ok")


rows = list(csv.DictReader(open(a.meta), delimiter="\t"))
changed = {r: collections.Counter() for r in RANKS}
status = collections.Counter()
blank_before = {r: 0 for r in RANKS}; blank_after = {r: 0 for r in RANKS}
with open(a.out, "w") as o:
    o.write("\t".join(["id", "taxid"] + RANKS + ["cluster"] + [f"{r}_taxid" for r in RANKS] + ["label", "taxid_status"]) + "\n")
    for r in rows:
        try:
            t = int(r["taxid"] or 0)
        except ValueError:
            t = 0
        lin, st = lineage(t) if t else ({}, "no_taxid")
        status[st] += 1
        new = {}
        for k in RANKS:
            if st in ("ok", "merged"):
                new[k] = lin.get(k, (0, ""))
            else:
                new[k] = (0, r[k])           # keep the BFD label, no taxid
            if not r[k]: blank_before[k] += 1
            if not new[k][1]: blank_after[k] += 1
            if st in ("ok", "merged") and r[k] != new[k][1]:
                changed[k][(r[k], new[k][1])] += 1
        o.write("\t".join([r["id"], str(t)] + [new[k][1] for k in RANKS] + [clusters.get(r["id"], "")] +
                          [str(new[k][0]) for k in RANKS] + [r["species"], st]) + "\n")

n = len(rows)
md = [f"# Taxonomy normalization\n", f"{n} genomes; NCBI taxdump: {a.taxdump_label or a.taxdump}.\n",
      "Status of the BFD NCBI taxids: " + ", ".join(f"{k} {v}" for k, v in status.most_common()) + ".\n",
      "| rank | labels changed | blank before | blank after | most common changes (BFD -> NCBI) |", "|---|---|---|---|---|"]
for k in RANKS:
    tot = sum(changed[k].values())
    top = "; ".join(f"{b or '(blank)'} -> {c or '(blank)'} x{v}" for (b, c), v in changed[k].most_common(4))
    md.append(f"| {k} | {tot} ({100*tot/n:.1f}%) | {blank_before[k]} | {blank_after[k]} | {top} |")
if a.report:
    open(a.report, "w").write("\n".join(md) + "\n")
print("\n".join(md), file=sys.stderr)
