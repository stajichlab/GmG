#!/usr/bin/env python3
"""D1g: does the visiting order of the collapse change what is found downstream?
usage: order_check.py SLUG SEED BASELINE_DIR RUN_DIR >> order_check.tsv
BASELINE_DIR is a run without --collapse; RUN_DIR a run with --collapse (and --collapse-seed SEED, or `noise` for the default order).
Columns: slug seed n reps k_clusters ari_vs_baseline study_r2 study_p country_r2 country_p
ari_vs_baseline: adjusted Rand index between the baseline clusters of all samples and the clusters of their representatives
(every collapsed sample takes the cluster of its representative)."""
import csv
import json
import sys
from collections import Counter

import numpy as np


def clusters(path):
    with open(path) as f:
        r = csv.reader(f, delimiter="\t")
        next(r)
        return {row[0]: row[1] for row in r}


def ari(a, b):
    ka, kb = sorted(set(a)), sorted(set(b))
    ia, ib = {x: i for i, x in enumerate(ka)}, {x: i for i, x in enumerate(kb)}
    m = np.zeros((len(ka), len(kb)))
    for x, y in zip(a, b):
        m[ia[x], ib[y]] += 1
    c2 = lambda v: v * (v - 1) / 2
    s_ij = c2(m).sum()
    s_a, s_b = c2(m.sum(1)).sum(), c2(m.sum(0)).sum()
    n2 = c2(len(a))
    exp = s_a * s_b / n2 if n2 else 0.0
    mx = (s_a + s_b) / 2
    return 1.0 if mx == exp else (s_ij - exp) / (mx - exp)


def tests(path):
    out = {}
    try:
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                out[row["label_column"]] = (row["r2"], row["p"])
    except FileNotFoundError:
        pass
    return out


def main():
    slug, seed, base, run = sys.argv[1:5]
    cb = clusters(f"{base}/clusters.tsv")
    cr = clusters(f"{run}/clusters.tsv")
    with open(f"{run}/collapse.tsv") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    ids = [r["id"] for r in rows if r["id"] in cb]
    proj = [cr.get(next(r["representative"] for r in rows if r["id"] == i), "NA") for i in ids]
    reps = sum(r["is_representative"] == "true" for r in rows)
    k = len(set(cr.values()) - {"NA", "-", ""})
    t = tests(f"{run}/group_tests.tsv")
    g = lambda c, i: t.get(c, ("", ""))[i]
    print("\t".join(str(x) for x in [slug, seed, len(rows), reps, k, f"{ari([cb[i] for i in ids], proj):.3f}",
                                      g("study", 0), g("study", 1), g("country", 0), g("country", 1)]), flush=True)


if __name__ == "__main__":
    main()
