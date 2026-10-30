#!/usr/bin/env python3
"""Stability of the cluster ids between index releases (#47, #48).

usage: cluster_stability.py OUT.md LABEL=meta_full.tsv [LABEL=meta_full.tsv ...] [--excluded FILE]
The `cluster` column of a release's meta_full.tsv is the id of the cluster centroid (95% ANI). For every consecutive pair of releases, on
the genomes both contain: the share whose centroid id is unchanged, the number of clusters, the adjusted Rand index of the two
partitions, and why genomes moved (their old centroid left the reference, or was excluded, or the genome joined a different centroid).
"""
import collections, csv, sys


def load(p):
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    return {r["id"]: r["cluster"] for r in rows if r.get("cluster")}


def comb2(n):
    return n * (n - 1) // 2


def ari(a, b):
    cont = collections.Counter((a[i], b[i]) for i in a)
    ra, rb = collections.Counter(a.values()), collections.Counter(b[i] for i in a)
    sum_c = sum(comb2(v) for v in cont.values())
    sum_a = sum(comb2(v) for v in ra.values())
    sum_b = sum(comb2(v) for v in rb.values())
    n2 = comb2(len(a))
    exp = sum_a * sum_b / n2 if n2 else 0
    mx = (sum_a + sum_b) / 2
    return (sum_c - exp) / (mx - exp) if mx != exp else 1.0


def main():
    args = sys.argv[1:]
    out = args.pop(0)
    excl = set()
    if "--excluded" in args:
        i = args.index("--excluded")
        excl = {l.strip() for l in open(args[i + 1]) if l.strip()}
        del args[i:i + 2]
    rel = [(x.split("=", 1)[0], load(x.split("=", 1)[1])) for x in args]
    md = ["# Cluster-id stability between releases\n",
          "Cluster id = id of the 95% ANI centroid (`cluster` column of `meta_full.tsv`). Compared on the genomes both releases contain.\n",
          "| from | to | genomes in both | same centroid id | clusters (from / to) | adjusted Rand index | centroid left the reference | centroid excluded | joined another centroid |",
          "|---|---|---|---|---|---|---|---|---|"]
    for (la, a), (lb, b) in zip(rel, rel[1:]):
        common = {i: a[i] for i in a if i in b}
        same = sum(1 for i in common if a[i] == b[i])
        moved = [i for i in common if a[i] != b[i]]
        gone = sum(1 for i in moved if a[i] not in b)
        ex = sum(1 for i in moved if a[i] not in b and a[i] in excl)
        other = len(moved) - gone
        ca = len({a[i] for i in common})
        cb = len({b[i] for i in common})
        md.append(f"| {la} | {lb} | {len(common)} | {same} ({100 * same / len(common):.2f}%) | {ca} / {cb} | {ari(common, {i: b[i] for i in common}):.4f} | {gone} | {ex} | {other} |")
    md.append("")
    md.append("`centroid left the reference`: the old centroid is not in the new release (excluded, or dropped by the genome filter); its cluster is re-seeded by "
              "the next genome in the BFD priority order. `joined another centroid`: the centroid exists in both releases but the genome was assigned elsewhere "
              "(greedy clustering is order dependent: a new genome that outranks a centroid can take members).\n")
    open(out, "w").write("\n".join(md) + "\n")
    print("\n".join(md))


main()
