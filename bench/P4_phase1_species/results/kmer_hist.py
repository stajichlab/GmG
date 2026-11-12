#!/usr/bin/env python3
"""P4: the k-mer abundance histogram of some sketches (whole sketch, scaled 200): counts of hashes at abundance 1..16 and 17+.
usage: kmer_hist.py SKETCH_DIR IDS_FILE   (run on HPCC; needs numpy)"""
import struct
import sys

import numpy as np


def abunds(p):
    b = open(p, "rb").read()
    o = 8 + 4 + 4 + 8
    nl = struct.unpack_from("<I", b, o)[0]
    o += 4 + nl
    n = struct.unpack_from("<Q", b, o)[0]
    o += 8 + 8 * n
    assert b[o] == 1
    return np.frombuffer(b, dtype="<u4", count=n, offset=o + 1)


d, ids = sys.argv[1], [l.strip() for l in open(sys.argv[2]) if l.strip()]
print("id".ljust(12), " ".join(f"{k:>5d}" for k in range(1, 17)), "  17+   total")
for i in ids:
    a = abunds(f"{d}/{i}.gsk")
    h = np.bincount(np.minimum(a, 17), minlength=18)[1:]
    print(i.ljust(12), " ".join(f"{x:5d}" for x in h[:16]), f"{h[16]:5d} {len(a):7d}")
