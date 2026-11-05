#!/usr/bin/env python3
"""09_bd_eval.py OUT_DIR PATRISTIC.tsv : score k-mer structure against Bd lineages and the SNP/gene tree distances."""
import sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score as ari
d, pat = sys.argv[1:3]
A = pd.read_csv(f"{d}/ani_est.tsv", sep="\t"); A.index = A.columns
man = pd.read_csv("manifest.tsv", sep="\t").set_index("id")
ids = list(A.index); lab = man.organism.reindex(ids); note = man.notes.reindex(ids)
asg = lab[lab != "Bd_unassigned"].index
print(f"{len(ids)} samples; {len(asg)} with lineage")
# within/between lineage ANI
L = sorted(lab[asg].unique()); M = pd.DataFrame(index=L, columns=L, dtype=float)
for a in L:
    for b in L:
        ia, ib = lab.index[lab == a], lab.index[lab == b]
        v = A.loc[ia, ib].values
        if a == b: v = v[np.triu_indices(len(ia), 1)]
        M.loc[a, b] = np.median(v) if v.size else np.nan
print("median ANI between lineages:\n", M.round(4).to_string())
D = 1 - A.loc[asg, asg].values; np.fill_diagonal(D, 0); D = (D + D.T) / 2
Z = linkage(squareform(D, checks=False), "average")
for k in (2, 3, 4, 6, 8):
    print(f"ARI(UPGMA k={k}, lineage) = {ari(lab[asg].values, fcluster(Z, k, 'maxclust')):.3f}")
h = pd.read_csv(f"{d}/pca_hash.tsv", sep="\t").set_index("id")
for k in (3, 4, 8):
    km = KMeans(k, n_init=10, random_state=0).fit_predict(h.loc[asg, ["PC1", "PC2", "PC3"]])
    print(f"ARI(hash-PCA kmeans k={k}, lineage) = {ari(lab[asg].values, km):.3f}")
# nearest-neighbour lineage agreement (leave-one-out 1-NN)
Aa = A.loc[asg, asg].values.copy(); np.fill_diagonal(Aa, -1)
nn = np.array(asg)[Aa.argmax(1)]; print(f"1-NN lineage agreement = {(lab[nn].values == lab[asg].values).mean():.3f}")
# genus-only samples
go = note.index[note == "bd_genus_only"]
print(f"genus-only samples in set: {len(go)}; their median nearest-neighbour ANI = {np.median([A.loc[s].drop(s).max() for s in go]):.4f} vs others {np.median([A.loc[s].drop(s).max() for s in ids if s not in set(go)]):.4f}")
print("genus-only lineage mix:", lab[go].value_counts().to_dict())
# tree distance concordance
P = pd.read_csv(pat, sep="\t", index_col=0)
sh = [i for i in ids if i in P.index and i in P.columns]
iu = np.triu_indices(len(sh), 1)
pv = P.loc[sh, sh].values[iu]; av = (1 - A.loc[sh, sh].values)[iu]
ok = ~np.isnan(pv)
print(f"{len(sh)} samples in patristic matrix; Spearman(patristic, 1-ANI) = {spearmanr(pv[ok], av[ok]).correlation:.3f}")
cl = pv[ok] < np.quantile(pv[ok], 0.5)
print(f"  closest half of pairs: {spearmanr(pv[ok][cl], av[ok][cl]).correlation:.3f}")
