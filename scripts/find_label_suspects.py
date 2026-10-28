#!/usr/bin/env python3
"""Find genomes whose metadata lineage conflicts with their near-identical neighbours.

usage: find_label_suspects.py BENCH_LOO.tsv OUT_PREFIX [MIN_ANI=0.95] [MIN_CONS=0.8]
Uses the leave-one-out bench (every genome vs all others): a genome is a suspect at a rank if the vote of its
closest hits (top hits within 10% of the best shared count) names a *different* taxon with consistency >= MIN_CONS
and the best hit's query-containment ANI >= MIN_ANI. At that identity two genomes cannot be in different
phyla/classes/orders, so a conflict at those ranks is a label error (or contamination) in one of them.
Also lists "hub" references: genomes that are the top hit of >= 3 conflicting queries (the likely culprit when many
genomes of one group all disagree with one reference).
Writes OUT_PREFIX.suspects.tsv, OUT_PREFIX.hubs.tsv, OUT_PREFIX.exclude.txt (suspect queries + hubs) and prints counts.
"""
import sys
import pandas as pd

bench, out = sys.argv[1], sys.argv[2]
min_ani = float(sys.argv[3]) if len(sys.argv) > 3 else 0.95
min_cons = float(sys.argv[4]) if len(sys.argv) > 4 else 0.8
d = pd.read_csv(bench, sep="\t", dtype=str, keep_default_na=False)
for c in ["containment_query", "top_shared"] + [f"{r}_cons" for r in ["phylum", "class", "order", "family"]]:
    d[c] = pd.to_numeric(d[c])
d["ani"] = d["containment_query"] ** (1 / 21)
rows = []
for r in ["phylum", "class", "order"]:
    x = d[(d[f"{r}_true"] != "") & (d[f"{r}_vote"] != "") & (d[f"{r}_vote"] != d[f"{r}_true"]) & (d[f"{r}_cons"] >= min_cons) & (d["ani"] >= min_ani)]
    for _, row in x.iterrows():
        rows.append(dict(id=row["id"], rank=r, own=row[f"{r}_true"], neighbours=row[f"{r}_vote"], top_id=row["top_id"], ani=round(row["ani"], 4),
                         genus=row["genus_true"], species=row["species_true"]))
s = pd.DataFrame(rows)
print(f"{len(d)} genomes in the leave-one-out run; ANI >= {min_ani}, vote consistency >= {min_cons}")
if s.empty:
    sys.exit(0)
for r in ["phylum", "class", "order"]:
    print(f"  conflicts at {r}: {(s['rank']==r).sum()}")
s.to_csv(out + ".suspects.tsv", sep="\t", index=False)
cls = s[s["rank"].isin(["phylum", "class"])]
hubs = cls.groupby("top_id").agg(n=("id", "size"), neighbours_say=("neighbours", lambda v: v.value_counts().index[0]),
                                 their_label=("own", lambda v: v.value_counts().index[0])).query("n >= 3").sort_values("n", ascending=False)
hubs.to_csv(out + ".hubs.tsv", sep="\t")
print(f"  hub references (top hit of >=3 phylum/class-conflicting queries): {len(hubs)}")
excl = sorted(set(cls["id"]) | set(hubs.index))
open(out + ".exclude.txt", "w").write("\n".join(excl) + "\n")
print(f"  exclude list: {len(excl)} genomes")
print(hubs.head(8).to_string())
