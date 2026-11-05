#!/usr/bin/env python3
"""11_snp_pca.py : SNP PCA (from the existing alignments) vs k-mer hash PCA, same samples.

SNP PCA: biallelic sites from the alignment (ACGT only, minor allele in >= MINOR samples, <= 20% missing among the
compared samples), minor-allele dosage 0/1, missing -> column mean, centred, randomized PCA.
Comparison on shared samples: R^2 of each k-mer PC regressed on the SNP PCs (and back), principal-angle cosines of the
top-3/top-5 subspaces, Procrustes disparity (top 3), Spearman of pairwise distances in the top-5 PC space.
"""
import sys
import numpy as np, pandas as pd
from scipy.linalg import subspace_angles
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

G = "/bigdata/stajichlab/shared/projects/Onygenales/Popgen/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/strain_tree/unpruned/"
BD_ALN = "/bigdata/stajichlab/shared/projects/Chytrid/Bd_popgen/Phylogeny/BdStrains.347_taxa.JGI_1086.cds.fasaln"
W = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/"
MINOR = 5
LUT = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i; LUT[b + 32] = i


def read_aln(path, want=None):
    seqs, name, buf = {}, None, []
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if name is not None and (want is None or name in want): seqs[name] = LUT[np.frombuffer(b"".join(buf), dtype=np.uint8)]
                name, buf = line[1:].split()[0].decode(), []
            else:
                buf.append(line.strip())
        if name is not None and (want is None or name in want): seqs[name] = LUT[np.frombuffer(b"".join(buf), dtype=np.uint8)]
    return seqs


def snp_pca(seqs, names, ncomp=10):
    X = np.vstack([seqs[n] for n in names])               # n x L, 0-3 or 255
    n = len(names)
    cnt = np.stack([(X == b).sum(0) for b in range(4)])   # 4 x L
    order = np.argsort(-cnt, axis=0)
    major = order[0]; minor = order[1]
    cmaj = np.take_along_axis(cnt, major[None], 0)[0]; cmin = np.take_along_axis(cnt, minor[None], 0)[0]
    called = cnt.sum(0)
    keep = (cmin >= MINOR) & (called >= 0.8 * n)
    Xk = X[:, keep]; mj = major[keep]; mn = minor[keep]
    D = np.full(Xk.shape, np.nan, dtype=np.float32)
    D[Xk == mj] = 0; D[Xk == mn] = 1                      # third alleles / N / gaps stay missing
    cm = np.nanmean(D, axis=0); idx = np.where(np.isnan(D)); D[idx] = cm[idx[1]]
    D -= D.mean(0)
    p = PCA(n_components=ncomp, svd_solver="randomized", random_state=0)
    return p.fit_transform(D), p.explained_variance_ratio_, int(keep.sum())


def r2(y, X):
    X1 = np.column_stack([np.ones(len(X)), X]); b, *_ = np.linalg.lstsq(X1, y, rcond=None)
    return 1 - ((y - X1 @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()


def compare(tag, K, S, labels, ax_row, evr, nsites):
    """K: n x 5 k-mer PCs; S: n x 10 SNP PCs (same sample order)."""
    print(f"\n=== {tag}: n={len(K)} shared samples, {nsites} SNP sites; SNP PC1-3 var = {np.round(evr[:3], 3).tolist()}")
    print("  R2 of k-mer PC_j from SNP PC1-5 / PC1-10:", [f"PC{j+1}: {r2(K[:, j], S[:, :5]):.2f}/{r2(K[:, j], S[:, :10]):.2f}" for j in range(5)])
    print("  R2 of SNP  PC_j from k-mer PC1-5       :", [f"PC{j+1}: {r2(S[:, j], K[:, :5]):.2f}" for j in range(5)])
    for m in (3, 5):
        c = np.cos(subspace_angles(K[:, :m], S[:, :m])); print(f"  principal-angle cosines top-{m}: {np.round(np.sort(c)[::-1], 3).tolist()}")
    Kz = (K[:, :3] - K[:, :3].mean(0)) / K[:, :3].std(0); Sz = (S[:, :3] - S[:, :3].mean(0)) / S[:, :3].std(0)
    print(f"  Procrustes disparity (top 3, 0=identical, 1=unrelated): {procrustes(Kz, Sz)[2]:.3f}")
    print(f"  Spearman of pairwise distances in top-5 PC space: {spearmanr(pdist(K[:, :5]), pdist(S[:, :5])).correlation:.3f}")
    pal = ["#d95f02", "#1b9e77", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#666666", "#a6761d"]
    cats = sorted(set(labels)); col = {c: pal[i % len(pal)] for i, c in enumerate(cats)}
    for a, M, nm in ((ax_row[0], K, "k-mer hash PCA"), (ax_row[1], S, "SNP PCA")):
        for c in cats:
            m = np.array(labels) == c; a.scatter(M[m, 0], M[m, 1], s=11, c=col[c], label=c, alpha=.85)
        a.set_title(f"{tag}: {nm} (n={len(K)})", fontsize=9); a.set_xlabel("PC1"); a.set_ylabel("PC2")
    ax_row[0].legend(fontsize=6, loc="best")


fig, ax = plt.subplots(4, 2, figsize=(10, 15))
# ---- Coccidioides
for r, (tag, kdir, mfa, kk) in enumerate([("Cocci all", "final_s200", "Cocciall_v1.All.SNP.mfa", 0),
                                         ("C. immitis", "final_s200_immitis", "Cocciall_v1.Cimmitis.SNP.mfa", 2),
                                         ("C. posadasii", "final_s200_posadasii", "Cocciall_v1.Cposadasii.SNP.mfa", 3)]):
    h = pd.read_csv(f"{W}popgen_cocci/{kdir}/pca_hash.tsv", sep="\t").set_index("id")
    g = pd.read_csv(f"{W}popgen_cocci/{kdir}/groups.tsv", sep="\t").set_index("id")
    seqs = read_aln(G + mfa, set(h.index))
    names = [i for i in h.index if i in seqs]
    S, evr, ns = snp_pca(seqs, names)
    K = h.loc[names, ["PC1", "PC2", "PC3", "PC4", "PC5"]].values
    labels = list(g.loc[names, "organism_label"].str.replace("Coccidioides ", "C. ")) if kk == 0 else \
        [f"k-mer cluster {c}" for c in KMeans(kk, n_init=10, random_state=0).fit_predict(K[:, :2])]
    compare(tag, K, S, labels, ax[r], evr, ns)
    pd.DataFrame(S[:, :5], index=names).to_csv(f"{W}popgen_cocci/snp_pca_{kdir}.tsv", sep="\t")
# ---- Bd (gene-CDS SNPs from the 347-taxon alignment)
man = pd.read_csv(f"{W}popgen_bd/manifest.tsv", sep="\t").set_index("id")
tm = pd.read_csv("/rhome/jstajich/projects/TE_fungi_copy/results/bd/relatedness/bd_tree_map.tsv", sep="\t").dropna(subset=["tree_tip"]).set_index("strain")
h = pd.read_csv(f"{W}popgen_bd/out_s200/pca_hash.tsv", sep="\t").set_index("id")
ids = [i for i in h.index if i in tm.index]
seqs = read_aln(BD_ALN, set(tm.loc[ids, "tree_tip"]))
ids = [i for i in ids if tm.tree_tip[i] in seqs]
S, evr, ns = snp_pca(seqs, [tm.tree_tip[i] for i in ids])
K = h.loc[ids, ["PC1", "PC2", "PC3", "PC4", "PC5"]].values
coarse = {"Bd_L2": "L2-4", "Bd_L3": "L2-4", "Bd_L4": "L2-4", "Bd_L5": "L5/7", "Bd_L7": "L5/7", "Bd_L1": "L1", "Bd_L6": "L6", "Bd_L8": "L8"}
labels = [coarse.get(man.organism[i], "unassigned") for i in ids]
compare("Bd (gene CDS SNPs)", K, S, labels, ax[3], evr, ns)
pd.DataFrame(S[:, :5], index=ids).to_csv(f"{W}popgen_bd/snp_pca.tsv", sep="\t")
plt.tight_layout(); plt.savefig(f"{W}popgen_cocci/kmer_vs_snp_pca.png", dpi=100)
