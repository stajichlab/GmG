#!/usr/bin/env python3
"""All-pairs sketch ANI for a set of GmG .gsk sketches (same estimator family as gmg: containment^(1/k)).

For genomes a, b: shared = |A n B| after downsampling to --scaled; ani = (shared / min(|A|, |B|))^(1/k)
(max-containment ANI, so a smaller or less complete assembly is not penalised). Pairs with no shared hash are omitted.

usage: genome_ani.py --scaled 10000 --threads 16 -o pairs.tsv sketch_dir/*.gsk   (or --list FILE)
output: a b n_a n_b shared ani
"""
import argparse
import math
import os
import struct
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

MAGIC = b"GMGSKT01"


def read_gsk(path):
    """Return (name, k, scaled, sorted uint64 hashes)."""
    with open(path, "rb") as f:
        if f.read(8) != MAGIC:
            raise ValueError(f"{path}: not a .gsk file")
        k, scaled = struct.unpack("<II", f.read(8))
        f.read(8)  # seed
        (nl,) = struct.unpack("<I", f.read(4))
        name = f.read(nl).decode()
        (n,) = struct.unpack("<Q", f.read(8))
        h = np.frombuffer(f.read(8 * n), dtype="<u8").astype(np.uint64)
    return name, k, scaled, np.sort(h)


def downsample(h, scaled):
    return h[h <= np.uint64((2**64 - 1) // scaled)]


def count_shared(a, b):
    """|a n b| for sorted unique arrays."""
    if len(a) > len(b):
        a, b = b, a
    if len(a) == 0:
        return 0
    i = np.searchsorted(b, a)
    i[i == len(b)] = len(b) - 1
    return int((b[i] == a).sum())


_H = []
_K = 21


def _init(hs, k):
    global _H, _K
    _H, _K = hs, k


def _row(i):
    out = []
    a = _H[i]
    for j in range(i + 1, len(_H)):
        s = count_shared(a, _H[j])
        if s:
            m = min(len(a), len(_H[j]))
            out.append((i, j, s, (s / m) ** (1.0 / _K)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sketches", nargs="*")
    ap.add_argument("--list")
    ap.add_argument("--scaled", type=int, default=10000)
    ap.add_argument("--threads", type=int, default=os.cpu_count())
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    paths = list(a.sketches)
    if a.list:
        paths += [l.strip() for l in open(a.list) if l.strip()]
    names, hs, k = [], [], None
    for p in paths:
        name, kk, sc, h = read_gsk(p)
        if sc > a.scaled:
            sys.exit(f"{p}: sketch scaled {sc} is coarser than --scaled {a.scaled}")
        k = k or kk
        if kk != k:
            sys.exit(f"{p}: k={kk}, expected {k}")
        names.append(os.path.basename(p)[: -len(".gsk")] if p.endswith(".gsk") else name)
        hs.append(downsample(h, a.scaled))
    print(f"{len(hs)} sketches, k={k}, scaled={a.scaled}, mean {np.mean([len(h) for h in hs]):.0f} hashes", file=sys.stderr)
    with ProcessPoolExecutor(a.threads, initializer=_init, initargs=(hs, k)) as ex, open(a.out, "w") as w:
        w.write("a\tb\tn_a\tn_b\tshared\tani\n")
        for rows in ex.map(_row, range(len(hs)), chunksize=8):
            for i, j, s, ani in rows:
                w.write(f"{names[i]}\t{names[j]}\t{len(hs[i])}\t{len(hs[j])}\t{s}\t{ani:.5f}\n")


if __name__ == "__main__":
    main()
