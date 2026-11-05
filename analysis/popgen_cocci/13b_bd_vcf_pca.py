#!/usr/bin/env python3
"""13b_bd_vcf_pca.py : genome-wide SNP PCA from the Bd VCF vs the k-mer hash PCA (same samples).

Needs 13a_bd_vcf_dump.sh output (vcf_dosage.tsv, vcf_samples.txt). Panel ids (SRR / Sample_ / strain names) are mapped
to VCF sample names through the TE_fungi_copy name files; unmatched samples are reported and left out."""
import re, sys, os
import numpy as np, pandas as pd
from scipy.linalg import subspace_angles
from scipy.spatial import procrustes
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

O = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_bd/"
T = "/rhome/jstajich/projects/TE_fungi_copy/"
PD = "/bigdata/stajichlab/shared/projects/Population_Genomics/B_dendrobatidis/"
norm = lambda s: re.sub(r"^bd", "", re.sub(r"[^a-z0-9]", "", str(s).lower()))

vs = [l.strip() for l in open(O + "vcf_samples.txt")]
vmap = {}
for v in vs: vmap.setdefault(norm(v), []).append(v)
man = pd.read_csv(O + "manifest.tsv", sep="\t").set_index("id")
meta = pd.read_csv(T + "results/bd/bd_metadata.tsv", sep="\t").set_index("strain")
tm = pd.read_csv(T + "results/bd/relatedness/bd_tree_map.tsv", sep="\t").set_index("strain")
sc = pd.read_csv(PD + "samples.csv"); srr2strain = dict(zip(sc.FileBase, sc.Strain))
key = {}
kf = "/bigdata/stajichlab/shared/projects/Chytrid/Bd_popgen/asm/samples.all.dat"
if os.path.exists(kf):
    for line in open(kf):
        p = line.rstrip("\n").split("\t") if "\t" in line else line.split()
        if len(p) >= 2: key[p[0]] = p[1]; key.setdefault(p[1], p[0])
def cands(p):
    out = [p, srr2strain.get(p), key.get(p)]
    if p in tm.index: out.append(tm.tree_tip.get(p))
    if p in meta.index: out += [meta.sheet_strain.get(p), meta.isolate.get(p)]
    return [c for c in out if isinstance(c, str) and c]
match, how = {}, {}
for p in man.index:
    for c in cands(p):
        h = vmap.get(norm(c))
        if h and len(h) == 1:
            match[p] = h[0]; how[p] = c; break
print(f"panel samples {len(man)}; matched to VCF {len(match)}; VCF samples {len(vs)}; VCF names claimed {len(set(match.values()))}")
dup = pd.Series(list(match.values())).value_counts(); dup = dup[dup > 1]
if len(dup): print("VCF names matched by >1 panel id (kept first):", dup.to_dict())
seen = set(); match = {p: v for p, v in match.items() if not (v in seen or seen.add(v))}

G = pd.read_csv(O + "vcf_dosage.tsv", sep="\t", header=None, dtype=np.int8).values   # sites x vcf samples
vi = {v: i for i, v in enumerate(vs)}
ids = [p for p in match]
X = G[:, [vi[match[p]] for p in ids]].astype(np.float32)                          # sites x n
X[X < 0] = np.nan
n = X.shape[1]; called = np.sum(~np.isnan(X), 1); af = np.nansum(X, 1) / (2 * np.maximum(called, 1))
maf = np.minimum(af, 1 - af)
keep = (called >= 0.8 * n) & (maf >= 0.02)
X = X[keep].T                                                                       # n x sites
cm = np.nanmean(X, 0); idx = np.where(np.isnan(X)); X[idx] = cm[idx[1]]
X -= X.mean(0)
print(f"SNPs used: {X.shape[1]} of {G.shape[0]} (MAF>=0.02, >=80% called among {n} matched samples)")
pca = PCA(n_components=10, svd_solver="randomized", random_state=0); S = pca.fit_transform(X)
print("SNP PCA variance PC1-5:", np.round(pca.explained_variance_ratio_[:5], 3).tolist())
pd.DataFrame(S[:, :5], index=ids).to_csv(O + "snp_pca_vcf.tsv", sep="\t")

h = pd.read_csv(O + os.environ.get("KDIR", "out_s200") + "/pca_hash.tsv", sep="\t").set_index("id")
sh = [i for i in ids if i in h.index]; ix = [ids.index(i) for i in sh]
K = h.loc[sh, ["PC1", "PC2", "PC3", "PC4", "PC5"]].values; S = S[ix]
def r2(y, M):
    M1 = np.column_stack([np.ones(len(M)), M]); b, *_ = np.linalg.lstsq(M1, y, rcond=None)
    return 1 - ((y - M1 @ b) ** 2).sum() / ((y - y.mean()) ** 2).sum()
print(f"\n=== Bd genome-wide SNPs vs k-mer hash PCA: n={len(sh)} shared samples")
print("  R2 of k-mer PC_j from SNP PC1-5 / PC1-10:", [f"PC{j+1}: {r2(K[:, j], S[:, :5]):.2f}/{r2(K[:, j], S[:, :10]):.2f}" for j in range(5)])
print("  R2 of SNP  PC_j from k-mer PC1-5       :", [f"PC{j+1}: {r2(S[:, j], K[:, :5]):.2f}" for j in range(5)])
for m in (3, 5): print(f"  principal-angle cosines top-{m}:", np.round(np.sort(np.cos(subspace_angles(K[:, :m], S[:, :m])))[::-1], 3).tolist())
Kz = (K[:, :3] - K[:, :3].mean(0)) / K[:, :3].std(0); Sz = (S[:, :3] - S[:, :3].mean(0)) / S[:, :3].std(0)
print(f"  Procrustes disparity (top 3): {procrustes(Kz, Sz)[2]:.3f}")
print(f"  Spearman of pairwise distances, top-5 PC space: {spearmanr(pdist(K[:, :5]), pdist(S[:, :5])).correlation:.3f}")
# compare against the earlier gene-CDS SNP PCA on common samples
g = pd.read_csv(O + "snp_pca.tsv", sep="\t", index_col=0)
c = [i for i in sh if i in g.index]
if len(c) > 20:
    A = g.loc[c].values[:, :5]; B = S[[sh.index(i) for i in c]][:, :5]; Kc = K[[sh.index(i) for i in c]]
    print(f"  (gene-CDS SNP PCA vs genome-wide VCF PCA on {len(c)} samples: top-3 cosines "
          f"{np.round(np.sort(np.cos(subspace_angles(A[:, :3], B[:, :3])))[::-1], 3).tolist()}; "
          f"dist Spearman {spearmanr(pdist(A), pdist(B)).correlation:.3f}; k-mer vs gene-CDS {spearmanr(pdist(Kc), pdist(A)).correlation:.3f})")
coarse = {"Bd_L2": "L2-4", "Bd_L3": "L2-4", "Bd_L4": "L2-4", "Bd_L5": "L5/7", "Bd_L7": "L5/7", "Bd_L1": "L1", "Bd_L6": "L6", "Bd_L8": "L8"}
lab = np.array([coarse.get(man.organism[i], "unassigned") for i in sh])
pal = {"L1": "#d95f02", "L2-4": "#1b9e77", "L5/7": "#7570b3", "L6": "#e7298a", "L8": "#66a61e", "unassigned": "#e6ab02"}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
for a, M, t in ((ax[0], K, "k-mer hash PCA"), (ax[1], S, "genome-wide SNP PCA (VCF)")):
    for c_ in pal:
        m = lab == c_; a.scatter(M[m, 0], M[m, 1], s=11, c=pal[c_], label=c_, alpha=.85)
    a.set_title(f"Bd {t} (n={len(sh)})", fontsize=9); a.set_xlabel("PC1"); a.set_ylabel("PC2")
ax[0].legend(fontsize=7); plt.tight_layout(); plt.savefig(O + os.environ.get("KDIR", "out_s200") + "_vs_vcf_pca.png", dpi=100)
