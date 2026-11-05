#!/usr/bin/env python3
"""30_bd_shoulder_het.py : are the Bd samples flagged `abundance_shoulder` heterozygous diploids or mixtures?

Two tests. (1) Heterozygosity from the Bd VCF (fraction of heterozygous calls among called sites per strain): are the flagged
samples more heterozygous than the rest? (2) Where does the shoulder sit? A heterozygous peak is at half the main coverage;
a second organism would be anywhere below it. Needs the popsketch run (popsketch_runs/bd_full), the dosage matrix and name
files made by 13a_bd_vcf_dump.sh and used by 13b_bd_vcf_pca.py."""
import os, re, struct
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu

P = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/"
T = "/rhome/jstajich/projects/TE_fungi_copy/"
PD = "/bigdata/stajichlab/shared/projects/Population_Genomics/B_dendrobatidis/"
norm = lambda s: re.sub(r"^bd", "", re.sub(r"[^a-z0-9]", "", str(s).lower()))
q = pd.read_csv(P + "popsketch_runs/bd_full/qc.tsv", sep="\t").set_index("id")
q["flags"] = q["flags"].fillna("")
shoulder = set(q.index[q["flags"].str.contains("abundance_shoulder")])
print(f"Bd samples: {len(q)}; flagged abundance_shoulder: {len(shoulder)}")

# ---- (1) heterozygosity from the VCF
vs = [l.strip() for l in open(P + "popgen_bd/vcf_samples.txt")]
vmap = {}
for v in vs:
    vmap.setdefault(norm(v), []).append(v)
man = pd.read_csv(P + "popgen_bd/manifest.tsv", sep="\t").set_index("id")
meta = pd.read_csv(T + "results/bd/bd_metadata.tsv", sep="\t").set_index("strain")
tm = pd.read_csv(T + "results/bd/relatedness/bd_tree_map.tsv", sep="\t").set_index("strain")
sc = pd.read_csv(PD + "samples.csv"); srr2strain = dict(zip(sc.FileBase, sc.Strain))
match = {}
for p in man.index:
    cands = [p, srr2strain.get(p)]
    if p in tm.index: cands.append(tm.tree_tip.get(p))
    if p in meta.index: cands += [meta.sheet_strain.get(p), meta.isolate.get(p)]
    for c in cands:
        if isinstance(c, str) and c and len(vmap.get(norm(c), [])) == 1:
            match[p] = vmap[norm(c)][0]; break
G = pd.read_csv(P + "popgen_bd/vcf_dosage.tsv", sep="\t", header=None, dtype=np.int8).values
vi = {v: i for i, v in enumerate(vs)}
het = {}
for p, v in match.items():
    col = G[:, vi[v]]
    called = col >= 0
    het[p] = (col[called] == 1).mean() if called.any() else np.nan
h = pd.Series(het)
print(f"mapped to the VCF: {len(h)} of {len(q)} samples; flagged and mapped: {len(shoulder & set(h.index))}")
fl = h[h.index.isin(shoulder)].dropna(); rest = h[~h.index.isin(shoulder)].dropna()
print(f"heterozygous fraction of called sites: flagged median {fl.median():.4f} (n={len(fl)}), others median {rest.median():.4f} (n={len(rest)})")
if len(fl) >= 3:
    u = mannwhitneyu(fl, rest, alternative="greater")
    print(f"  Mann-Whitney (flagged > others): p = {u.pvalue:.3g}; AUC = {u.statistic / (len(fl) * len(rest)):.3f}")
print("  flagged samples' het fraction:", {k: round(v, 4) for k, v in fl.items()})
print("  top of the rest:", {k: round(v, 4) for k, v in rest.sort_values(ascending=False).head(6).items()})

# ---- (2) where the shoulder sits relative to c/2
def abunds(path):
    with open(path, "rb") as f:
        f.read(8); k, scaled, seed, nl = struct.unpack("<IIQI", f.read(20)); f.read(nl)
        (n,) = struct.unpack("<Q", f.read(8)); f.seek(8 * n, 1); f.read(1)
        return np.frombuffer(f.read(4 * n), dtype="<u4")
rows = []
for sid in sorted(shoulder) + list(q.index[(q["flags"] == "") & (q.status == "included")][:40]):
    p = P + "popsketch_runs/bd_full/sketches/" + sid + ".gsk"
    c = q.kmer_coverage.get(sid)
    if not os.path.exists(p) or not np.isfinite(c) or c < 6: continue
    a = abunds(p); a = a[a >= q.abundance_cutoff[sid]]
    lo, hi = 0.1 * c, 0.6 * c
    sh = a[(a >= lo) & (a < hi)]
    if len(sh) < 50: continue
    # mass of the shoulder within +-20% of c/2 versus the rest of the shoulder window
    near_half = ((sh >= 0.4 * c) & (sh <= 0.6 * c)).mean()
    rows.append(dict(id=sid, flagged=sid in shoulder, cov=c, shoulder_hashes=len(sh), frac_near_half=near_half, mode=float(np.bincount(sh).argmax()) / c))
d = pd.DataFrame(rows).set_index("id")
print("\nshoulder position (fraction of the shoulder window's hashes at 0.4c-0.6c; mode of the shoulder in units of c):")
print(d.groupby("flagged")[["cov", "shoulder_hashes", "frac_near_half", "mode"]].median().round(3).to_string())
print(d[d.flagged].round(3).to_string())
