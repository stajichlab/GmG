#!/usr/bin/env python3
"""D1: how many representatives are left if samples above an ANI threshold are collapsed?

Reads popsketch's ani.tsv (square matrix, ANI as a fraction) and prints, per threshold, the number of
 - leader clusters: walk samples in file order; each unassigned sample becomes a representative and takes every
   unassigned sample with ANI >= t to it (no chaining; every member is within t of its representative)
 - linked groups: connected components of the graph "ANI >= t" (chains allowed; an upper bound on collapse)
usage: collapse_counts.py SLUG ani.tsv  >> results.tsv   (columns: slug n threshold leaders components largest_leader)
"""
import sys
import numpy as np

THRESHOLDS = [0.99, 0.995, 0.999, 0.9995, 0.9999]


def load(path):
    with open(path) as f:
        f.readline()
        rows = [line.rstrip("\n").split("\t")[1:] for line in f]
    return np.array(rows, dtype=float)


def leaders(a, t):
    n = a.shape[0]
    free = np.ones(n, bool)
    sizes = []
    for i in range(n):
        if not free[i]:
            continue
        take = free & (a[i] >= t)
        take[i] = True
        sizes.append(int(take.sum()))
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


def main():
    slug, path = sys.argv[1], sys.argv[2]
    a = load(path)
    for t in THRESHOLDS:
        sz = leaders(a, t)
        print(f"{slug}\t{a.shape[0]}\t{t}\t{len(sz)}\t{components(a, t)}\t{max(sz)}", flush=True)


if __name__ == "__main__":
    main()
