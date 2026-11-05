#!/usr/bin/env python3
"""03_plots.py OUT_DIR -> OUT_DIR/overview.png (ANI histogram, PCoA, DAPC, genome size vs coverage)."""
import sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, pandas as pd
d = sys.argv[1]
col = {"Coccidioides immitis": "#d95f02", "Coccidioides posadasii": "#1b9e77", "Fusarium solani": "#7570b3"}
g = pd.read_csv(f"{d}/groups.tsv", sep="\t"); pc = pd.read_csv(f"{d}/pcoa.tsv", sep="\t").merge(g, on="id")
q = pd.read_csv(f"{d}/qc.tsv", sep="\t"); h = pd.read_csv(f"{d}/ani_hist.tsv", sep="\t")
var = [float(x) for x in open(f"{d}/pcoa_var.txt").read().split()]
fig, ax = plt.subplots(2, 2, figsize=(11, 9))
ax[0, 0].bar(h.ani_lo, h.pairs, width=0.0023, color="#555"); ax[0, 0].set_yscale("log")
ax[0, 0].set_xlabel("pairwise ANI (containment^(1/21))"); ax[0, 0].set_ylabel("pairs"); ax[0, 0].set_title("pairwise ANI")
for lab, c in col.items():
    s = pc[pc.organism_label == lab]; ax[0, 1].scatter(s.PC1, s.PC2, s=12, c=c, label=lab, alpha=.8)
ax[0, 1].legend(fontsize=7); ax[0, 1].set_xlabel(f"PC1 {var[0]:.0%}"); ax[0, 1].set_ylabel(f"PC2 {var[1]:.0%}"); ax[0, 1].set_title("PCoA (1-ANI)")
s = pc[pc.organism_label != "Fusarium solani"]
for lab, c in list(col.items())[:2]:
    t = s[s.organism_label == lab]; ax[1, 0].scatter(t.PC2, t.PC3, s=12, c=c, alpha=.8)
ax[1, 0].set_xlabel(f"PC2 {var[1]:.0%}"); ax[1, 0].set_ylabel(f"PC3 {var[2]:.0%}"); ax[1, 0].set_title("PCoA, Coccidioides only")
for lab, c in col.items():
    t = q[q.organism == lab]; ax[1, 1].scatter(t.median_kmer_cov, t.est_genome_mb, s=12, c=c, alpha=.8)
ax[1, 1].set_xscale("log"); ax[1, 1].set_xlabel("median k-mer coverage"); ax[1, 1].set_ylabel("est. size (Mb) = kept hashes x scaled")
ax[1, 1].axhspan(25, 30, color="#bbb", alpha=.3); ax[1, 1].set_title("QC: size vs coverage")
plt.tight_layout(); plt.savefig(f"{d}/overview.png", dpi=110)
