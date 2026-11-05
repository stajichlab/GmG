#!/usr/bin/env python3
"""14_bd_kmer_pc1.py : what does Bd k-mer hash PC1 track (coverage, genome size, ploidy flags, lineage)?"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
O = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_bd/"
h = pd.read_csv(O + "out_s200/pca_hash.tsv", sep="\t").set_index("id")
q = pd.read_csv(O + "out_s200/qc.tsv", sep="\t").set_index("id")
a = pd.read_csv("/rhome/jstajich/projects/TE_fungi_copy/resources/panels/bd_species_audit.tsv", sep="\t").set_index("strain")
m = pd.read_csv(O + "manifest.tsv", sep="\t").set_index("id")
d = h.join(q[["est_genome_mb", "median_kmer_cov", "kept_hashes"]]).join(a[["gb", "fail_reason", "sc_frac_present"]]).join(m[["organism"]])
for c in ("est_genome_mb", "median_kmer_cov", "kept_hashes", "gb", "sc_frac_present"):
    for pc in ("PC1", "PC2"):
        print(f"Spearman({pc}, {c}) = {spearmanr(d[pc], d[c], nan_policy='omit').correlation:.3f}")
d["gtl"] = (d.fail_reason == "genome_too_large")
print(d.groupby("gtl")[["PC1", "PC2", "est_genome_mb"]].median().round(1).to_string())
# natural groups along PC1
d["pc1_bin"] = pd.cut(d.PC1, [-300, -100, -50, 0, 100])
print(d.groupby("pc1_bin", observed=True).agg(n=("PC1", "size"), cov=("median_kmer_cov", "median"), mb=("est_genome_mb", "median"),
      gtl=("gtl", "mean"), L1=("organism", lambda x: (x == "Bd_L1").sum()), L8=("organism", lambda x: (x == "Bd_L8").sum()),
      L24=("organism", lambda x: x.isin(["Bd_L2", "Bd_L3", "Bd_L4"]).sum())).round(2).to_string())
