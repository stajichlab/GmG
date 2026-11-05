#!/usr/bin/env python3
"""25_hist_dump.py SK.gsk ... : abundance histogram (first 70 bins, then coarse) for each sketch, to see where a cutoff goes wrong."""
import struct, sys, os
import numpy as np
for p in sys.argv[1:]:
    with open(p, "rb") as f:
        f.read(8); k, scaled, seed, nl = struct.unpack("<IIQI", f.read(20)); f.read(nl)
        (n,) = struct.unpack("<Q", f.read(8)); f.seek(8 * n, 1); f.read(1); a = np.frombuffer(f.read(4 * n), dtype="<u4")
    h = np.bincount(np.minimum(a, 3000))
    print(f"{os.path.basename(p)[:-4]}: n={n}, max={a.max()}, mean={a.mean():.1f}, median={np.median(a):.0f}")
    print("  a=1..40:", h[1:41].tolist())
    print("  a=41..120 step 4 (sums):", [int(h[i:i + 4].sum()) for i in range(41, 121, 4)])
