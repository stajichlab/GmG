#!/usr/bin/env python3
"""D1e: collapse with a per-pair threshold from each sample's own split-half distance.
usage: per_pair.py SLUG ani.tsv split_half.tsv >> per_pair.tsv
A sample j joins representative i when 1 - ANI(i, j) <= tau(i, j). Rules:
  mean  tau = (d_i + d_j) / 2         (the typical noise of the comparison)
  max   tau = max(d_i, d_j)           (looser: the noisier sample sets the limit)
  species-median  tau = median(d)     (D1c, for comparison)
d_i is the sample's split-half distance; a sample without halves gets the species median.
Columns: slug n median rule reps largest_cluster"""
import sys

import numpy as np

from collapse_counts import load


def read_d(path, ids):
    d = {}
    with open(path) as f:
        f.readline()
        for line in f:
            i, v = line.rstrip("\n").split("\t")
            d[i] = float(v)
    med = float(np.median(list(d.values())))
    return np.array([d.get(i, med) for i in ids]), med


def reps(dist, tau):
    n = dist.shape[0]
    free = np.ones(n, bool)
    sizes = []
    for i in range(n):
        if not free[i]:
            continue
        take = free & (dist[i] <= tau[i])
        take[i] = True
        sizes.append(int(take.sum()))
        free &= ~take
    return sizes


def main():
    slug, ani, sh = sys.argv[1:4]
    a = load(ani)
    with open(ani) as f:
        f.readline()
        ids = [line.split("\t", 1)[0] for line in f]
    d, med = read_d(sh, ids)
    dist = 1.0 - a
    rules = {
        "mean": (d[:, None] + d[None, :]) / 2,
        "max": np.maximum(d[:, None], d[None, :]),
        "species-median": np.full(dist.shape, med),
    }
    for name, tau in rules.items():
        sz = reps(dist, tau)
        print(f"{slug}\t{a.shape[0]}\t{med:.5f}\t{name}\t{len(sz)}\t{max(sz)}", flush=True)


if __name__ == "__main__":
    main()
