#!/usr/bin/env python3
"""Clade-aware calibration strategies (issue #58), compared on held-out genera.

usage: clade_calibration.py BENCH_DIR EXCLUDE_IDS.txt OUT_PREFIX [N_SPLITS=10] [TARGET=0.95]
  BENCH_DIR has bench_{loo,los,log}_s10000.tsv from `gmg bench`; GMG_META selects the reference metadata (as analyze_bench.py).
Pool = `known` (the default thresholds' pool). Strategies, each fitted on half the genera and scored on the other half:
  global    one (min_ani, min_consistency) per rank at TARGET (what the shipped thresholds do)
  margin    the same fitted at TARGET + 0.02
  weighted  one fit per rank where every phylum carries equal total weight
  phylum    a fit per PREDICTED phylum (the query's own phylum vote) when it has >= MIN_N voted queries, else global
Scored per TRUE phylum group (Ascomycota / Basidiomycota / other), because the problem is under-delivery outside Ascomycota.
"""
import importlib.util, os, sys
import numpy as np
import pandas as pd

bench, excl, out = sys.argv[1:4]
n_splits = int(sys.argv[4]) if len(sys.argv) > 4 else 10
target = float(sys.argv[5]) if len(sys.argv) > 5 else 0.95
MIN_N = 400
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ab", os.path.join(here, "analyze_bench.py"))
ab = importlib.util.module_from_spec(spec)
argv, sys.argv = sys.argv, ["x"]
spec.loader.exec_module(ab)
sys.argv = argv
RANKS = ab.RANKS
POOL = {"species": ["loo"], "genus": ["loo"], "family": ["los"], "order": ["los"], "class": ["los"], "phylum": ["los"]}
CONS = (0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

d = ab.load([os.path.join(bench, f"bench_{m}_s10000.tsv") for m in ("loo", "los", "log")])
d = d[~d["id"].isin(set(open(excl).read().split()))].reset_index(drop=True)
ab.NC, ab.NSP, ab.NGE = ab.counts()
d = ab.add_answerable(d)
d["group"] = np.where(d["genus_true"] != "", d["genus_true"], np.where(d["species_true"] != "", d["species_true"], d["id"]))
d["pgroup"] = np.where(d["phylum_true"].isin(["Ascomycota", "Basidiomycota"]), d["phylum_true"], "other")


def frame(dd, rank):
    x = ab.rank_frame(dd[dd["mode"].isin(POOL[rank])], rank)
    return x[x["voted"]]


def fit(x, tgt, w=None):
    """(ani, cons) maximising (weighted) coverage subject to (weighted) precision >= tgt; None if no fit."""
    if len(x) < 50:
        return None
    w = np.ones(len(x)) if w is None else np.asarray(w, float)
    ani, cons, ok = x["ani"].values, x["cons"].values, x["vote_ok"].values.astype(float)
    grid = np.unique(np.round(np.percentile(ani, np.arange(0, 100, 1.0)), 4))
    best, tot = None, w.sum()
    for c in CONS:
        mc = cons >= c
        for a in grid:
            m = mc & (ani >= a)
            if m.sum() < 30:
                break
            prec = (w[m] * ok[m]).sum() / w[m].sum()
            cov = w[m].sum() / tot
            if prec >= tgt and (best is None or cov > best[2]):
                best = (float(a), c, cov)
    return best


def clade_weights(x):
    """Each phylum gets equal total weight (counted on its true-phylum group)."""
    g = x["phylum_true"]
    cnt = g.map(g.value_counts())
    return (1.0 / cnt).values * len(x) / g.nunique()


def apply(x, params, per_phylum=None):
    """Boolean mask of queries that would be called, given global params or per-predicted-phylum params."""
    if per_phylum is None:
        if params is None:
            return np.zeros(len(x), bool)
        return ((x["ani"] >= params[0]) & (x["cons"] >= params[1])).values
    m = np.zeros(len(x), bool)
    pv = x["phylum_vote"].values
    for ph, p in per_phylum.items():
        sel = pv == ph if ph != "*" else ~np.isin(pv, [k for k in per_phylum if k != "*"])
        if p is not None:
            m |= sel & (x["ani"].values >= p[0]) & (x["cons"].values >= p[1])
    return m


rows = []
groups = d["group"].unique()
for s in range(n_splits):
    rng = np.random.default_rng(2000 + s)
    cal_g = set(rng.permutation(groups)[: len(groups) // 2])
    cal, tst = d[d["group"].isin(cal_g)], d[~d["group"].isin(cal_g)]
    for rank in RANKS:
        xc, xt = frame(cal, rank), frame(tst, rank)
        if len(xc) < 100 or len(xt) < 100:
            continue
        g = fit(xc, target)
        strat = {
            "global": (apply(xt, g), None),
            "margin": (apply(xt, fit(xc, target + 0.02)), None),
            "weighted": (apply(xt, fit(xc, target, clade_weights(xc))), None),
        }
        vc = xc["phylum_vote"].value_counts()
        pp = {ph: fit(xc[xc["phylum_vote"] == ph], target) for ph in vc[vc >= MIN_N].index}
        pp["*"] = g
        strat["phylum"] = (apply(xt, None, pp), None)
        ok = xt["vote_ok"].values
        for name, (m, _) in strat.items():
            for grp in ["all", "Ascomycota", "Basidiomycota", "other"]:
                sel = np.ones(len(xt), bool) if grp == "all" else (xt["pgroup"].values == grp)
                n_sel = (m & sel).sum()
                rows.append(dict(split=s, rank=rank, strategy=name, group=grp, n=int(sel.sum()), called=int(n_sel),
                                 correct=int(ok[m & sel].sum()), cov=n_sel / max(sel.sum(), 1)))
    print("split", s, "done", file=sys.stderr)

r = pd.DataFrame(rows)
g = r.groupby(["rank", "strategy", "group"]).agg(called=("called", "sum"), correct=("correct", "sum"), n=("n", "sum")).reset_index()
g["precision"] = g["correct"] / g["called"].clip(lower=1)
g["coverage"] = g["called"] / g["n"].clip(lower=1)
g.to_csv(out + ".tsv", sep="\t", index=False, float_format="%.4f")

md = [f"# Clade-aware calibration strategies (target {target})\n",
      f"{len(d)} bench queries, {n_splits} random genus splits; fit on half the genera, score on the other half, pooled over splits. Pool `known` (leave-one-out for genus/species, leave-species-out above). "
      f"`phylum` fits a threshold per predicted phylum when it has ≥ {MIN_N} voted queries.\n"]
for grp in ["all", "Ascomycota", "Basidiomycota", "other"]:
    md.append(f"\n## True phylum group: {grp}\n")
    md.append("| rank | " + " | ".join(f"{s} prec (cov)" for s in ["global", "margin", "weighted", "phylum"]) + " |\n|---|---|---|---|---|")
    for rank in RANKS:
        cells = []
        for s in ["global", "margin", "weighted", "phylum"]:
            q = g[(g["rank"] == rank) & (g.strategy == s) & (g.group == grp)]
            cells.append("-" if q.empty else f"{q.precision.iloc[0]:.3f} ({q.coverage.iloc[0]:.2f})")
        md.append(f"| {rank} | " + " | ".join(cells) + " |")
open(out + ".md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
