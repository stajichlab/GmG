#!/usr/bin/env python3
"""02a_refhashes.py SAMPLES_CSV SKETCH_DIR OUT.npy [GENUS] -> union of reference hashes of one genus (ids -> SKETCH_DIR/id.gsk)."""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from importlib import import_module
rd = import_module("02_popdist").read_gsk
csv_path, skdir, out = sys.argv[1:4]; genus = sys.argv[4] if len(sys.argv) > 4 else "Coccidioides"
ids = [r["ASMID"] for r in csv.DictReader(open(csv_path)) if r["GENUS"] == genus]
hs, miss, sc = [], [], set()
for i in ids:
    p = os.path.join(skdir, i + ".gsk")
    if not os.path.exists(p): miss.append(i); continue
    k, s, _, h, _ = rd(p); sc.add((k, s)); hs.append(h)
u = np.unique(np.concatenate(hs)); np.save(out, u)
print(f"{len(hs)} sketches ({len(miss)} missing) k/scaled={sc} union={len(u)} -> {out}")
