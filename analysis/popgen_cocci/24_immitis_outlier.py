#!/usr/bin/env python3
"""24_immitis_outlier.py : who dominates hash PC1 in the immitis rerun, and is the 19-sample cluster still there without it?"""
import numpy as np, pandas as pd
P = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci/"
h = pd.read_csv(P + "final_s200_immitis_rc/pca_hash.tsv", sep="\t", index_col=0)
q = pd.read_csv(P + "final_s200_immitis_rc/qc.tsv", sep="\t").set_index("id")
nb = pd.read_csv(P + "final_s200_immitis_rc/neighbors.tsv", sep="\t").set_index("id")
z = (h.PC1 - h.PC1.median()).abs().sort_values(ascending=False)
print("largest |PC1| samples:"); print(h.loc[z.index[:4], ["PC1", "PC2", "PC3"]].round(1).join(q[["median_kmer_cov", "est_genome_mb", "min_abund"]]).join(nb[["nn1_ani"]]).round(3).to_string())
old = pd.read_csv(P + "final_s200_immitis/pca_hash.tsv", sep="\t", index_col=0)
x = old.index[old.PC1 > 40]
print(f"\nearlier 19-sample cluster: {len(x)} samples")
print("their PC1/PC2/PC3 in the corrected-cutoff PCA (median vs the rest):")
rest = h.index.difference(x)
print(pd.DataFrame({"cluster": h.loc[x, ["PC1", "PC2", "PC3"]].median(), "rest": h.loc[rest, ["PC1", "PC2", "PC3"]].median(),
                    "rest_sd": h.loc[rest, ["PC1", "PC2", "PC3"]].std()}).round(1).to_string())
