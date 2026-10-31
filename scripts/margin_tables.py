#!/usr/bin/env python3
"""Before/after tables for the genus/family safety margin (#58) from two validate_run.py .tsv files per dataset.

usage: margin_tables.py OUT.md NAME=BEFORE.tsv,AFTER.tsv [NAME=...]
Per dataset, per group (all, Ascomycota, Basidiomycota, other phyla pooled) and per rank: n with truth, called share, precision
with its 95% Wilson interval, before -> after. Counts are rebuilt from the table (called = round(called x n), correct = round(precision x called)).
"""
import sys
import math
import csv

RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (c - h, c + h)


def counts(path):
    out = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        g = r["group"]
        if g != "all" and not g.startswith("phylum="):
            continue
        n = int(float(r["n"]))
        called = round(float(r["called"]) * n)
        correct = round(float(r["precision"] or 0) * called) if called else 0
        if g.startswith("phylum="):
            ph = g.split("=", 1)[1]
            g = ph if ph in ("Ascomycota", "Basidiomycota") else "other phyla"
        o = out.setdefault((g, r["rank"]), [0, 0, 0])
        o[0] += n; o[1] += called; o[2] += correct
    return out


def cell(o):
    n, c, k = o
    if n == 0:
        return "-"
    if c == 0:
        return f"0 called (n={n})"
    lo, hi = wilson(k, c)
    return f"{c / n:.3f} / {k / c:.3f} [{lo:.3f}, {hi:.3f}] (n={n})"


md = []
for spec in sys.argv[2:]:
    name, files = spec.split("=", 1)
    b, a = (counts(f) for f in files.split(","))
    md.append(f"\n## {name}\n\nCells: called share / precision [95% Wilson CI] (n with truth). The CI is on the calls made.\n")
    for g in ["all", "Ascomycota", "Basidiomycota", "other phyla"]:
        md.append(f"\n**{g}**\n\n| rank | before | after | precision change | called change |\n|---|---|---|---|---|")
        for r in RANKS:
            ob, oa = b.get((g, r), [0, 0, 0]), a.get((g, r), [0, 0, 0])
            if ob[0] == 0 and oa[0] == 0:
                continue
            pb = ob[2] / ob[1] if ob[1] else float("nan")
            pa = oa[2] / oa[1] if oa[1] else float("nan")
            cb = ob[1] / ob[0] if ob[0] else float("nan")
            ca = oa[1] / oa[0] if oa[0] else float("nan")
            md.append(f"| {r} | {cell(ob)} | {cell(oa)} | {pa - pb:+.3f} | {ca - cb:+.3f} |")
open(sys.argv[1], "w").write("\n".join(md) + "\n")
