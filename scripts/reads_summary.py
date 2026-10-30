#!/usr/bin/env python3
"""One row per run from several `gmg validate` outputs (validation.tsv), for the reads validations (#42/#43/#45).

usage: reads_summary.py GROUP LABEL=validation.tsv [LABEL=validation.tsv ...]
GROUP is a `group` value of validation.tsv (`all`, `genus_in_ref=True`, `frac=10`, ...). Prints a Markdown table: for each run (in the order given)
the number of inputs with truth and "called / precision" at phylum, class, order, family and genus (and species where the truth has it).
Cells are `called share / precision (calls)`.
"""
import csv, sys

group = sys.argv[1]
runs = [a.split("=", 1) for a in sys.argv[2:]]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
out = ["| run | n | " + " | ".join(RANKS) + " |", "|---|---|" + "---|" * len(RANKS)]
for label, path in runs:
    d = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["group"] == group: d[r["rank"]] = r
    if not d:
        out.append(f"| {label} | 0 |" + " - |" * len(RANKS)); continue
    n = max(int(r["n"]) for r in d.values())
    cells = []
    for k in RANKS:
        r = d.get(k)
        if not r: cells.append("-"); continue
        called = float(r["called"]); calls = round(called * int(r["n"]))
        prec = r["precision"]
        cells.append(f"{called:.2f} / {float(prec):.3f} ({calls})" if prec not in ("", "nan") else f"{called:.2f} / - (0)")
    out.append(f"| {label} | {n} | " + " | ".join(cells) + " |")
print("\n".join(out))
