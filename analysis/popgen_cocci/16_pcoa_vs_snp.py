#!/usr/bin/env python3
"""PCoA(1-ANI) vs hash-presence PCA, each compared with SNP PCA, same shared samples."""
import os
import numpy as np, pandas as pd
from scipy.linalg import subspace_angles
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr

SUF = os.environ.get("SUF", "")  # e.g. _rc for the runs that use the Rust cutoffs
C = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci/"
B = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_bd/"
sets = [("Cocci all", C + "final_s200" + SUF + "/", C + "snp_pca_final_s200.tsv"),
        ("C. immitis", C + "final_s200_immitis" + SUF + "/", C + "snp_pca_final_s200_immitis.tsv"),
        ("C. posadasii", C + "final_s200_posadasii" + SUF + "/", C + "snp_pca_final_s200_posadasii.tsv"),
        ("Bd (VCF, cleaned)", B + "out_s200_qc" + SUF + "/", B + "snp_pca_vcf.tsv")]


def metrics(K, S):
    c3 = np.sort(np.cos(subspace_angles(K[:, :3], S[:, :3])))[::-1]
    c5 = np.sort(np.cos(subspace_angles(K[:, :5], S[:, :5])))[::-1]
    z = lambda M: (M[:, :3] - M[:, :3].mean(0)) / M[:, :3].std(0)
    return (f"cos3 {np.round(c3, 3).tolist()} | cos5 min {c5.min():.3f} | Procrustes {procrustes(z(K), z(S))[2]:.3f} | "
            f"dist rho {spearmanr(pdist(K[:, :5]), pdist(S[:, :5])).correlation:.3f}")


for name, d, snp in sets:
    S = pd.read_csv(snp, sep="\t", index_col=0)
    P = pd.read_csv(d + "pcoa.tsv", sep="\t", index_col=0)
    H = pd.read_csv(d + "pca_hash.tsv", sep="\t", index_col=0)
    ids = [i for i in S.index if i in P.index and i in H.index]
    s = S.loc[ids].values[:, :5]
    print(f"\n== {name}: n={len(ids)}")
    print("  hash PCA :", metrics(H.loc[ids].values[:, :5], s))
    print("  PCoA     :", metrics(P.loc[ids].values[:, :5], s))
