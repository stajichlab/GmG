#!/usr/bin/env python3
"""Where would 'cut UPGMA at the largest gap in merge heights' cut, on each dataset? Top gaps, k they imply, ratio gap1/gap2."""
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

C = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci/"
B = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_bd/"
for name, d in [("Cocci all (2 species)", C + "final_s200/"), ("C. immitis only", C + "final_s200_immitis/"),
                ("C. posadasii only", C + "final_s200_posadasii/"), ("Bd cleaned", B + "out_s200_qc/")]:
    A = pd.read_csv(d + "ani_est.tsv", sep="\t"); D = 1 - A.values; np.fill_diagonal(D, 0); D = (D + D.T) / 2
    n = len(D); h = linkage(squareform(D, checks=False), "average")[:, 2]     # n-1 merge heights, ascending
    # cut between merge i and i+1 gives k = n - 1 - i clusters; gap = h[i+1]-h[i]; look at the top 6 cuts, k<=30
    gaps = [(h[i + 1] - h[i], n - 1 - i, h[i], h[i + 1]) for i in range(len(h) - 1) if n - 1 - i <= 30]
    gaps.sort(reverse=True)
    print(f"\n== {name}: n={n}, top merge height {h[-1]:.4f}")
    for g, k, lo, hi in gaps[:5]:
        print(f"   k={k:2d}  gap {g:.4f}  (merge heights {lo:.4f} -> {hi:.4f})")
    print(f"   ratio gap1/gap2 = {gaps[0][0] / gaps[1][0]:.2f}  -> rule picks k={gaps[0][1]}")
