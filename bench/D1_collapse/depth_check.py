#!/usr/bin/env python3
"""D1d: distance of a sample to its own lower-depth subsample, next to the split-half distance.
usage: depth_check.py SLUG ani.tsv >> depth_self.tsv
ani.tsv comes from a pooled run whose ids are <run>__<depth> (for example SRR1__s1g). For every run present at two depths,
prints the distance (1 - ANI) between the two sketches of the same run. Columns: slug run depth_a depth_b distance"""
import sys
from itertools import combinations

from collapse_counts import load


def main():
    slug, path = sys.argv[1:3]
    a = load(path)
    with open(path) as f:
        f.readline()
        ids = [line.split("\t", 1)[0] for line in f]
    by = {}
    for i, x in enumerate(ids):
        run, depth = x.rsplit("__", 1)
        by.setdefault(run, {})[depth] = i
    for run, d in sorted(by.items()):
        for da, db in combinations(sorted(d), 2):
            print(f"{slug}\t{run}\t{da}\t{db}\t{1.0 - a[d[da], d[db]]:.5f}", flush=True)


if __name__ == "__main__":
    main()
