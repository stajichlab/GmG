#!/usr/bin/env python3
"""D1: how many representatives are left if samples above an ANI threshold are collapsed?

Reads popsketch's ani.tsv (square matrix, ANI as a fraction) and prints, per threshold, the number of
 - leader clusters: walk samples in file order; each unassigned sample becomes a representative and takes every
   unassigned sample with ANI >= t to it (no chaining; every member is within t of its representative)
 - linked groups: connected components of the graph "ANI >= t" (chains allowed; an upper bound on collapse)
usage: collapse_counts.py SLUG ani.tsv  >> results.tsv   (columns: slug n threshold leaders components largest_leader)

With a noise floor (popsketch.json "noise_floor.value": the 90th percentile of split-half 1 - ANI), the thresholds are
1 - m * floor for m in FLOOR_MULTIPLES, and the representatives for one multiple are written to a file:
usage: collapse_counts.py SLUG ani.tsv --floor F --reps-at M --reps-out reps.txt  >> results_floor.tsv
(columns: slug n floor multiple threshold leaders components largest_leader)
"""
import argparse
import json
import sys
import numpy as np

FLOOR_MULTIPLES = [1, 2, 3, 5, 10]

THRESHOLDS = [0.99, 0.995, 0.999, 0.9995, 0.9999]


def load(path):
    with open(path) as f:
        f.readline()
        rows = [line.rstrip("\n").split("\t")[1:] for line in f]
    return np.array(rows, dtype=float)


def leaders(a, t, reps=None):
    n = a.shape[0]
    free = np.ones(n, bool)
    sizes = []
    for i in range(n):
        if not free[i]:
            continue
        take = free & (a[i] >= t)
        take[i] = True
        sizes.append(int(take.sum()))
        if reps is not None:
            reps.append(i)
        free &= ~take
    return sizes


def components(a, t):
    n = a.shape[0]
    seen = np.zeros(n, bool)
    c = 0
    for s in range(n):
        if seen[s]:
            continue
        c += 1
        stack = [s]
        seen[s] = True
        while stack:
            i = stack.pop()
            nb = np.where((a[i] >= t) & ~seen)[0]
            seen[nb] = True
            stack.extend(nb.tolist())
    return c


def ids_of(path):
    with open(path) as f:
        f.readline()
        return [line.split("\t", 1)[0] for line in f]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("slug")
    ap.add_argument("ani")
    ap.add_argument("--floor", type=float)
    ap.add_argument("--reps-at", type=float, help="floor multiple whose representatives are written")
    ap.add_argument("--reps-out")
    g = ap.parse_args()
    a = load(g.ani)
    if g.floor is None:
        for t in THRESHOLDS:
            sz = leaders(a, t)
            print(f"{g.slug}\t{a.shape[0]}\t{t}\t{len(sz)}\t{components(a, t)}\t{max(sz)}", flush=True)
        return
    ids = ids_of(g.ani)
    for m in FLOOR_MULTIPLES:
        t = 1.0 - m * g.floor
        reps = []
        sz = leaders(a, t, reps)
        print(f"{g.slug}\t{a.shape[0]}\t{g.floor:.5f}\t{m}\t{t:.5f}\t{len(sz)}\t{components(a, t)}\t{max(sz)}", flush=True)
        if g.reps_out and m == g.reps_at:
            with open(g.reps_out, "w") as f:
                f.write("\n".join(ids[i] for i in reps) + "\n")


if __name__ == "__main__":
    main()
