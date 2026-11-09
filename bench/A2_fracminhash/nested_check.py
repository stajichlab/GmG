#!/usr/bin/env python3
"""A2, check: is the sampling noise of the scaled-200 -> scaled-1000 step what the model in predict.py says?
usage: nested_check.py SLUG SKETCH_DIR [N_SAMPLES] [MIN_ABUND] >> nested_check.tsv
A scaled-1000 sketch is the scaled-200 sketch restricted to hashes <= 2^64/1000, so the two are nested and share their noise: for
a pair with C = ANI^k and M hashes at scaled 200, SD(C_1000 - C_200) = sqrt(C (1 - C) (4 / M)) (the variance of the extra
sampling of M/5 hashes, with the nested covariance subtracted). Samples are the first N (sorted ids); hashes with abundance
< MIN_ABUND are dropped in both (a fixed stand-in for popsketch's cutoff: the check is of the sampling, not of the cutoff).
Columns: slug n_pairs median_one_minus_ani rms_z frac_abs_z_gt_2 median_hashes_200 median_hashes_1000 bytes_ratio
z = (ANI_1000 - ANI_200) / predicted SD of that difference; the model holds when rms_z is near 1."""
import glob
import os
import struct
import sys
from itertools import combinations

import numpy as np

K = 21
FULL = 2.0**64


def read_gsk(path):
    b = open(path, "rb").read()
    o = 8 + 4 + 4 + 8
    nl = struct.unpack_from("<I", b, o)[0]
    o += 4 + nl
    n = struct.unpack_from("<Q", b, o)[0]
    o += 8
    hashes = np.frombuffer(b, dtype="<u8", count=n, offset=o)
    o += 8 * n
    assert b[o] == 1
    o += 1
    ab = np.frombuffer(b, dtype="<u4", count=n, offset=o)
    return hashes, ab


def containment(a, b):
    """C = |a & b| / |a|, ANI-style symmetric: geometric mean of the two containments"""
    shared = np.intersect1d(a, b, assume_unique=True).size
    return shared / max(len(a), 1), shared / max(len(b), 1), shared


def ani(c1, c2):
    c = math.sqrt(max(c1, 0) * max(c2, 0))
    return c ** (1 / K), c


import math  # noqa: E402


def main():
    slug, d = sys.argv[1:3]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 80
    min_ab = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    files = sorted(f for f in glob.glob(os.path.join(d, "*.gsk")) if "half" not in f)[:n]
    s200, s1000 = [], []
    for f in files:
        h, ab = read_gsk(f)
        h = np.sort(h[ab >= min_ab])
        s200.append(h)
        s1000.append(h[h <= np.uint64(FULL / 1000)])
    dist, zs = [], []
    for i, j in combinations(range(len(files)), 2):
        a2, c2 = ani(*containment(s200[i], s200[j])[:2])
        a1, _ = ani(*containment(s1000[i], s1000[j])[:2])
        m = (len(s200[i]) + len(s200[j])) / 2
        c = c2
        if c <= 0 or c >= 1:
            continue
        sd_c = math.sqrt(c * (1 - c) * 4 / m)
        sd_ani = (c2 ** (1 / K) / (K * c)) * sd_c
        dist.append(1 - a2)
        zs.append((a1 - a2) / sd_ani)
    zs = np.array(zs)
    h200 = np.median([len(x) for x in s200])
    h1000 = np.median([len(x) for x in s1000])
    print("\t".join(str(x) for x in [slug, len(zs), f"{np.median(dist):.5f}", f"{np.sqrt(np.mean(zs**2)):.2f}",
                                      f"{np.mean(np.abs(zs) > 2):.3f}", int(h200), int(h1000), f"{h200/h1000:.2f}"]), flush=True)


if __name__ == "__main__":
    main()
