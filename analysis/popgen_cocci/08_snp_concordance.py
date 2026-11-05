#!/usr/bin/env python3
"""08_snp_concordance.py : test sketch-based groups/PCA clusters against the existing SNP trees.

For each species: (1) Spearman correlation of SNP-tree patristic distance vs our 1-ANI over shared samples;
(2) for each hash-PCA k-means cluster: is it a clade of the SNP tree (exact), and the best-matching split (F1);
(3) ARI between the SNP tree cut into k clades (average linkage on patristic distance) and our clusters.
Species/label check uses the all-species tree."""
import sys, itertools
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from Bio import Phylo

BASE = "/bigdata/stajichlab/shared/projects/Onygenales/Popgen/2025_All_Cocci/Genotyping/all_C_immitis_ref_RS/strain_tree/unpruned/"
TREES = {"immitis": BASE + "Cocciall_v1.Cimmitis.SNP.mfa.treefile",
         "posadasii": BASE + "Cocciall_v1.Cposadasii.SNP.fasttree.tre",
         "all": BASE + "Cocciall_v1.All.SNP.mfa.treefile"}
OUT = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_cocci/"


def clade_sets(t):
    sets = {}
    def rec(c):
        if c.is_terminal():
            s = frozenset([c.name])
        else:
            s = frozenset().union(*[rec(x) for x in c.clades])
        sets[c] = s
        return s
    rec(t.root)
    return sets


def patristic(t, names):
    parent = {}
    for c in t.find_clades(order="level"):
        for x in c.clades:
            parent[x] = (c, x.branch_length or 0.0)
    term = {c.name: c for c in t.get_terminals()}
    anc = {}
    for n in names:
        d, node, acc = {}, term[n], 0.0
        d[node] = 0.0
        while node in parent:
            p, bl = parent[node]; acc += bl; d[p] = acc; node = p
        anc[n] = d
    N = len(names); M = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            a, b = anc[names[i]], anc[names[j]]
            M[i, j] = M[j, i] = min(a[x] + b[x] for x in (a.keys() & b.keys()))
    return M


def split_scores(sets, universe, X):
    """exact: X (or its complement) is a clade; also best F1 of any clade/complement vs X, within `universe`."""
    X = frozenset(X); comp = universe - X; best = (0, None); exact = False
    for c, s in sets.items():
        s = s & universe
        if not s or s == universe:
            continue
        for tgt in (X,):
            for cand in (s, universe - s):
                tp = len(cand & tgt)
                if tp == 0: continue
                f1 = 2 * tp / (len(cand) + len(tgt))
                if f1 > best[0]: best = (f1, len(cand))
                if cand == tgt: exact = True
    return exact, best


trees = {k: Phylo.read(v, "newick") for k, v in TREES.items()}
for k, t in trees.items():
    print(f"tree {k}: {len(t.get_terminals())} tips")

# ---- species / label check on the all-species tree
allt = trees["all"]; asets = clade_sets(allt); tips_all = {c.name for c in allt.get_terminals()}
lab = {}
for sp in ("immitis", "posadasii"):
    for i in pd.read_csv(f"{OUT}final_s200_{sp}/groups.tsv", sep="\t").id: lab[i] = sp
common = tips_all & set(lab)
print(f"\n== species (all-species SNP tree): {len(common)} of {len(lab)} sketch samples are in the tree")
for sp in ("immitis", "posadasii"):
    X = {i for i in common if lab[i] == sp}
    ex, (f1, n) = split_scores(asets, frozenset(common), X)
    print(f"{sp}: n={len(X)} exact_clade={ex} best_split_F1={f1:.4f}")
sw = {"485B-1_L_OLD_CPA0023": "posadasii (relabelled; manifest says immitis)", "B3476": "immitis (relabelled; manifest says posadasii)"}
for s, w in sw.items():
    print(f"{s}: in SNP tree={s in tips_all}; sketch label={w}")
# where do the swaps sit in the SNP tree: nearest SNP-tree neighbours
if tips_all >= set(sw):
    names = sorted(common); P = patristic(allt, names); idx = {n: i for i, n in enumerate(names)}
    for s in sw:
        if s in idx:
            o = np.argsort(P[idx[s]])[1:6]
            print(f"  {s} nearest in SNP tree: " + ", ".join(f"{names[j]}({lab[names[j]]},{P[idx[s], j]:.4f})" for j in o))
    for sp in ("immitis", "posadasii"):
        pass

# ---- within-species concordance
rows = []
K_POS = int(sys.argv[1]) if len(sys.argv) > 1 else 4
for sp, k in (("immitis", 2), ("posadasii", K_POS)):
    t = trees[sp]; sets = clade_sets(t); tips = {c.name for c in t.get_terminals()}
    d = f"{OUT}final_s200_{sp}/"
    h = pd.read_csv(d + "pca_hash.tsv", sep="\t").set_index("id")
    h["cl"] = KMeans(k, n_init=10, random_state=0).fit_predict(h[["PC1", "PC2"]])
    A = pd.read_csv(d + "ani_est.tsv", sep="\t"); A.index = A.columns
    shared = sorted(tips & set(h.index))
    print(f"\n== {sp}: {len(shared)} shared of {len(h)} sketch samples / {len(tips)} SNP-tree tips")
    P = patristic(t, shared)
    iu = np.triu_indices(len(shared), 1)
    ani = A.loc[shared, shared].values
    rho = spearmanr(P[iu], (1 - ani)[iu]).correlation
    print(f"Spearman(SNP patristic, 1-ANI) = {rho:.3f}  (all pairs)")
    # restrict to close pairs (within the same coarse region) to avoid long-branch domination
    close = P[iu] < np.quantile(P[iu], 0.5)
    print(f"Spearman, closest half of pairs = {spearmanr(P[iu][close], (1 - ani)[iu][close]).correlation:.3f}")
    cl = h.loc[shared, "cl"].values
    uni = frozenset(shared)
    for c in range(k):
        X = [s for s, x in zip(shared, cl) if x == c]
        ex, (f1, n) = split_scores(sets, uni, X)
        med = h.loc[X, ["PC1"]].median().iloc[0]
        print(f"cluster {c} (PC1~{med:.0f}): n={len(X)} exact_clade={ex} best_split_F1={f1:.3f}")
    D = squareform(P, checks=False)
    Z = linkage(D, "average")
    for kk in sorted({k, 2, 3, 5, 8}):
        print(f"  ARI(tree cut k={kk}, our {k} clusters) = {adjusted_rand_score(fcluster(Z, kk, 'maxclust'), cl):.3f}")
    # which SNP-tree-cut cluster each of our clusters sits in
    tc = fcluster(Z, k, "maxclust")
    print(pd.crosstab(pd.Series(cl, name="our"), pd.Series(tc, name=f"snp_tree_k{k}")).to_string())
    qc = pd.read_csv(d + "qc.tsv", sep="\t").set_index("id")
    ct = pd.crosstab(cl, tc); major = {c: ct.loc[c].idxmax() for c in ct.index}
    disc = [(s_, c_, t_) for s_, c_, t_ in zip(shared, cl, tc) if t_ != major[c_]]
    print(f"discordant (cluster's majority SNP clade differs): {len(disc)}/{len(shared)}")
    for s_, c_, t_ in disc: print(f"   {s_}: our cl {c_} -> snp clade {t_}; cov {qc.median_kmer_cov[s_]:.0f}x, {qc.est_genome_mb[s_]:.1f} Mb")
    pd.DataFrame({"id": shared, "our_cluster": cl, "snp_cut": tc}).to_csv(f"{OUT}snp_concordance_{sp}" + (f"_k{k}" if sp == "posadasii" else "") + ".tsv", sep="\t", index=False)
