#!/usr/bin/env python3
"""D1c: a collapse threshold per species, from that species' own split-half distances (all samples, not 20).

usage: species_threshold.py SLUG ani.tsv popsketch.json >> species_threshold.tsv
Candidate rules (distance = 1 - ANI):  half_med = 0.5 x median split-half distance, med = median, p90 = 90th percentile.
Columns: slug n floor_n med p90 nn_p10 nn_p50 nn_p90 reps_half_med reps_med reps_p90
 nn_*: quantiles of each sample's distance to its nearest neighbour (what the closest relative looks like);
 reps_*: leader-cluster representatives left at 1 - rule (counts from collapse_counts.leaders)."""
import json
import sys

import numpy as np

from collapse_counts import leaders, load


def main():
    slug, ani, js = sys.argv[1:4]
    a = load(ani)
    f = json.load(open(js))["noise_floor"]
    d = 1.0 - a
    np.fill_diagonal(d, np.inf)
    nn = d.min(axis=1)
    med, p90 = f["median"], f["p90"]
    out = [slug, a.shape[0], f["n"], f"{med:.5f}", f"{p90:.5f}"]
    out += [f"{np.quantile(nn, q):.5f}" for q in (0.1, 0.5, 0.9)]
    out += [len(leaders(a, 1.0 - r)) for r in (0.5 * med, med, p90)]
    print("\t".join(str(x) for x in out), flush=True)


if __name__ == "__main__":
    main()
