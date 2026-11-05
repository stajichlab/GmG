#!/usr/bin/env python3
"""10_bd_coarse.py OUT_DIR : merge Bd lineages the sketches cannot tell apart (2/3/4 and 5/7) and re-score."""
import sys, numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score as ari
A = pd.read_csv(f"{sys.argv[1]}/ani_est.tsv", sep="\t"); A.index = A.columns
m = pd.read_csv("manifest.tsv", sep="\t").set_index("id"); lab = m.organism.reindex(A.index)
coarse = {"Bd_L2": "A", "Bd_L3": "A", "Bd_L4": "A", "Bd_L5": "B", "Bd_L7": "B", "Bd_L1": "L1", "Bd_L6": "L6", "Bd_L8": "L8"}
asg = lab[lab != "Bd_unassigned"].index; cg = lab[asg].map(coarse)
D = 1 - A.loc[asg, asg].values; np.fill_diagonal(D, 0); D = (D + D.T) / 2
Z = linkage(squareform(D, checks=False), "average")
for k in (3, 4, 5, 6):
    print(f"ARI(UPGMA k={k}, coarse groups) = {ari(cg.values, fcluster(Z, k, 'maxclust')):.3f}")
c = fcluster(Z, 5, "maxclust"); print(pd.crosstab(pd.Series(cg.values, name="coarse"), pd.Series(c, name="UPGMA k=5")).to_string())
Aa = A.loc[asg, asg].values.copy(); np.fill_diagonal(Aa, -1); nn = np.array(asg)[Aa.argmax(1)]
print("1-NN coarse-group agreement =", round((cg.reindex(nn).values == cg.values).mean(), 3))
un = lab.index[lab == "Bd_unassigned"]; asgl = list(asg)
pred = [cg[asgl[int(np.argmax(A.loc[s, asgl].values))]] for s in un]
print("unassigned -> nearest coarse group:", pd.Series(pred).value_counts().to_dict())
