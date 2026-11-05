#!/usr/bin/env python3
"""18_cutoff_compare.py : gmg-core's auto_min_abund (Rust probe output) vs the valley-below-main-peak rule in 02_popdist.py.

usage: 18_cutoff_compare.py RUST.tsv PY_QC.tsv [id id ...]  (ids: samples to print individually)"""
import sys
import numpy as np, pandas as pd
r = pd.read_csv(sys.argv[1], sep="\t").set_index("id")
p = pd.read_csv(sys.argv[2], sep="\t").set_index("id")
d = r.join(p[["min_abund", "kept_hashes", "est_genome_mb", "median_kmer_cov"]], how="inner", rsuffix="_py")
d = d.rename(columns={"est_mb": "rust_mb", "est_genome_mb": "py_mb", "min_abund": "py_min_abund"})
med = d.loc[d.rust_min_abund == d.py_min_abund, "py_mb"].median()   # cohort size where both rules agree
d["rust_dev"] = (d.rust_mb - med).abs() / med; d["py_dev"] = (d.py_mb - med).abs() / med
print(f"{len(d)} samples; size where rules agree: median {med:.1f} Mb")
print(f"identical cutoff: {(d.rust_min_abund == d.py_min_abund).sum()}  | rust higher: {(d.rust_min_abund > d.py_min_abund).sum()}  | rust lower: {(d.rust_min_abund < d.py_min_abund).sum()}")
print("rust cutoff distribution:", d.rust_min_abund.value_counts().sort_index().to_dict())
print("python cutoff distribution:", d.py_min_abund.value_counts().sort_index().to_dict())
bad = d[(d.rust_dev > 0.5) | (d.py_dev > 0.5)]
print(f"\nsamples where either rule gives a size >50% from the cohort median: {len(bad)}")
diff = d[d.rust_min_abund != d.py_min_abund]
print(f"samples where the cutoffs differ: {len(diff)}; rust closer to cohort size in {(diff.rust_dev < diff.py_dev).sum()}, python closer in {(diff.py_dev < diff.rust_dev).sum()}")
cols = ["rust_min_abund", "py_min_abund", "rust_mb", "py_mb", "median_kmer_cov", "valley"]
if len(diff): print("\ndifferences (first 25, largest gap in size first):"); print(diff.assign(g=(diff.rust_mb - diff.py_mb).abs()).sort_values("g", ascending=False)[cols].head(25).round(1).to_string())
ids = [i for i in sys.argv[3:] if i in d.index]
if ids: print(f"\nnamed samples ({len(ids)} of {len(sys.argv[3:])} found):"); print(d.loc[ids, cols].round(1).to_string())
