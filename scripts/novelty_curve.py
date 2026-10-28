#!/usr/bin/env python3
"""How often does the tool name a taxon that is NOT in the reference, and what does a stricter cutoff cost? (BFD benchmark only)

usage: novelty_curve.py BENCH_DIR EXCLUDE_IDS OUT_PREFIX
  genus:   leave-one-out (genus has relatives) gives coverage/precision; leave-genus-out (no relative of the genus exists) gives the
           FALSE-CALL rate: any genus named there is wrong by construction.
  species: same with leave-one-out vs leave-species-out.
For a grid of minimum ANI cutoffs (vote consistency fixed at the deployed value). Choose the cutoff on BFD, then check it on an independent set.
"""
import sys, os
import numpy as np, pandas as pd

bench, excl, out = sys.argv[1:4]
ex = set(open(excl).read().split())
def load(m):
    d = pd.read_csv(os.path.join(bench, f"bench_{m}_s10000.tsv"), sep="\t", dtype=str, keep_default_na=False)
    d = d[~d["id"].isin(ex)].copy()
    d["ani"] = pd.to_numeric(d["containment_query"]) ** (1 / 21)
    for r in ("genus", "species"):
        d[f"{r}_cons"] = pd.to_numeric(d[f"{r}_cons"])
    return d
loo, los, log = load("loo"), load("los"), load("log")
rows = []
for rank, hard, cons, grid in (("genus", log, 0.5, [0.7765, 0.80, 0.85, 0.88, 0.90, 0.92, 0.95, 0.97]),
                               ("species", los, 0.9, [0.97, 0.98, 0.99, 0.9938, 0.997, 0.999])):
    base = loo[(loo[f"{rank}_true"] != "")]
    # answerable in leave-one-out: another genome carries the same taxon
    cnt = base[f"{rank}_true"].map(base[f"{rank}_true"].value_counts())
    base = base[cnt >= 2]
    h = hard[hard[f"{rank}_true"] != ""]
    for t in grid:
        sel = base[(base[f"{rank}_vote"] != "") & (base[f"{rank}_cons"] >= cons) & (base["ani"] >= t)]
        fc = h[(h[f"{rank}_vote"] != "") & (h[f"{rank}_cons"] >= cons) & (h["ani"] >= t)]
        rows.append(dict(rank=rank, min_ani=t, loo_coverage=len(sel) / len(base), loo_precision=(sel[f"{rank}_vote"] == sel[f"{rank}_true"]).mean(),
                         absent_taxon_false_call_rate=len(fc) / len(h), n_loo=len(base), n_absent=len(h)))
df = pd.DataFrame(rows)
df.to_csv(out + ".tsv", sep="\t", index=False, float_format="%.4f")
md = ["# False calls for taxa absent from the reference (BFD benchmark)\n",
      "`loo` = leave-one-out with relatives present: coverage and precision of the call. `absent` = the genus (or species) was left out of the reference, so any such call is wrong: the false-call rate is the share of those queries that still get a call. "
      f"Clean labels, {len(ex)} suspect genomes dropped.\n"]
for rank in ("genus", "species"):
    md.append(f"\n## {rank}\n\n| min ANI | leave-one-out coverage | leave-one-out precision | false-call rate when the {rank} is absent |\n|---|---|---|---|")
    for _, r in df[df["rank"] == rank].iterrows():
        md.append(f"| {r.min_ani:.4f} | {r.loo_coverage:.3f} | {r.loo_precision:.3f} | {r.absent_taxon_false_call_rate:.3f} |")
open(out + ".md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
