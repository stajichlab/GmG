#!/usr/bin/env python3
"""Cluster purity and cluster-level accuracy.

usage: analyze_clusters.py CLUSTERS.tsv META.tsv [BENCH_LOO.tsv]  > report.md
CLUSTERS.tsv from `gmg cluster` (id, cluster, is_centroid, ani_to_centroid, ...).
Reports cluster sizes, how well species labels agree with clusters (purity), how fragmented species are across
clusters, the most mixed clusters, and, with a leave-one-out bench, how often the best hit is in the query's own cluster
(vs. carrying the query's species label) for queries whose cluster has >= 2 members.
"""
import sys
import pandas as pd

cl = pd.read_csv(sys.argv[1], sep="\t", dtype=str, keep_default_na=False)
meta = pd.read_csv(sys.argv[2], sep="\t", dtype=str, keep_default_na=False)
d = cl.merge(meta[["id", "species", "genus", "class", "phylum"]].rename(columns=lambda c: c + "_m" if c != "id" else c), on="id")
d["species_m"] = d["species_m"].where(d["species_m"] != "", "(no species)")
sizes = d.groupby("cluster").size()
print(f"# cluster analysis\n\n{len(d)} genomes, {len(sizes)} clusters; singletons: {(sizes==1).sum()} ({(sizes==1).sum()/len(sizes):.0%} of clusters, {(sizes==1).sum()/len(d):.0%} of genomes); largest: {sizes.max()}\n")
q = sizes.quantile([.5, .9, .99]).tolist()
print(f"cluster size median / p90 / p99: {q[0]:.0f} / {q[1]:.0f} / {q[2]:.0f}\n")

# purity: share of genomes carrying their cluster's most common species label
maj = d.groupby("cluster")["species_m"].agg(lambda s: s.value_counts().iloc[0])
purity = maj.sum() / len(d)
print(f"**purity** (genomes whose species label is the most common label in their cluster): {purity:.3f}")
multi = d[d["cluster"].map(sizes) >= 2]
maj2 = multi.groupby("cluster")["species_m"].agg(lambda s: s.value_counts().iloc[0]).sum() / len(multi)
print(f"purity over clusters with >= 2 genomes: {maj2:.3f}\n")

# fragmentation: for each species with >= 5 genomes, share of its genomes in its largest cluster
sp_sizes = d.groupby("species_m").size()
big = sp_sizes[(sp_sizes >= 5) & (sp_sizes.index != "(no species)")].index
frag = d[d["species_m"].isin(big)].groupby("species_m")["cluster"].agg(lambda s: s.value_counts().iloc[0] / len(s))
print(f"species with >= 5 genomes: {len(big)}; median share of a species in its largest cluster: {frag.median():.2f}; "
      f"species split over >1 substantial cluster (largest holds < 80%): {(frag < 0.8).sum()}\n")

# mixed clusters
rows = []
for c, g in d.groupby("cluster"):
    vc = g["species_m"].value_counts()
    if len(vc) >= 2 and vc.iloc[1] >= 5:
        rows.append((len(g), c, "; ".join(f"{k} ({v})" for k, v in vc.head(4).items())))
rows.sort(reverse=True)
print(f"clusters in which a second species label has >= 5 genomes: {len(rows)}\n")
print("| size | centroid | labels |\n|---|---|---|")
for n, c, lab in rows[:15]:
    print(f"| {n} | {c} | {lab} |")

if len(sys.argv) > 3:
    b = pd.read_csv(sys.argv[3], sep="\t", dtype=str, keep_default_na=False)
    cmap = dict(zip(d["id"], d["cluster"]))
    b["c_true"] = b["id"].map(cmap); b["c_top"] = b["top_id"].map(cmap)
    b["csize"] = b["c_true"].map(sizes)
    ans = b[b["csize"] >= 2]
    sp_ok = (ans["species_vote"] == ans["species_true"])
    cl_ok = (ans["c_top"] == ans["c_true"])
    print(f"\n## leave-one-out, queries whose cluster has >= 2 members (n={len(ans)})\n")
    print(f"- best hit is in the query's own cluster: {cl_ok.mean():.3f}")
    print(f"- species label of the vote equals the query's label: {sp_ok.mean():.3f}")
    print(f"- of the queries with a wrong species label, the best hit is still in the query's cluster: {cl_ok[~sp_ok].mean():.3f} (n={int((~sp_ok).sum())})")

if "n_hashes" in cl.columns:
    cl["n_hashes"] = pd.to_numeric(cl["n_hashes"])
    dd = cl.merge(meta[["id", "species"]].rename(columns={"species": "sp"}), on="id")
    dd = dd[dd["sp"] != ""]
    med = dd.groupby("sp")["n_hashes"].transform("median"); cnt = dd.groupby("sp")["n_hashes"].transform("size")
    big = dd[(cnt >= 5) & (dd["n_hashes"] > 1.6 * med)].copy()
    big["x_median"] = (big["n_hashes"] / med[big.index]).round(2)
    print(f"\n## genome-size outliers\n\n{len(big)} genomes are > 1.6x the median size (sketch hashes) of their species (species with >= 5 genomes): possible mixed/contaminated assemblies, or true polyploid/hybrid genomes.\n")
    print("| id | species | x median | hashes |\n|---|---|---|---|")
    for _, r in big.sort_values("x_median", ascending=False).head(25).iterrows():
        print(f"| {r['id']} | {r['sp']} | {r['x_median']} | {r['n_hashes']} |")
    big[["id", "sp", "x_median", "n_hashes"]].to_csv("size_outliers.tsv", sep="\t", index=False)
