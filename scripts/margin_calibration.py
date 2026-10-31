#!/usr/bin/env python3
"""Safety margin at genus and family (issue #58): derive the cutoffs and score them on held-out genera.

usage: margin_calibration.py BENCH_DIR EXCLUDE_IDS.txt OUT_PREFIX [N_SPLITS=10] [TARGET=0.95] [MARGIN=0.02] --shipped-json known.json
  BENCH_DIR has bench_{loo,los,log}_s10000.tsv from `gmg bench`; GMG_META selects the reference metadata (as analyze_bench.py).
Pool `known` (genus/species leave-one-out, family and above leave-species-out), same protocol as clade_calibration.py:
  fit on half the genera, score on the other half, 10 random splits (seeds 2000+s), pooled.
Strategies per rank:
  shipped   the cutoffs in the --shipped-json file (data/thresholds/bfd_s10000.known.json before the margin)
  fit       a fit at TARGET on the calibration half (the clade_calibration `global`)
  margin    a fit at TARGET + MARGIN on the calibration half (the clade_calibration `margin`)
  new       the stricter of shipped and margin fit, per rank: the margin may only tighten a cutoff, never loosen one that was
            raised on purpose (genus ANI 0.90, #57)
Also writes the fits on the FULL pool (the numbers for the threshold file) at TARGET and TARGET + MARGIN.
"""
import argparse, importlib.util, json, os, sys
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("bench")
ap.add_argument("excl")
ap.add_argument("out")
ap.add_argument("n_splits", nargs="?", type=int, default=10)
ap.add_argument("target", nargs="?", type=float, default=0.95)
ap.add_argument("margin", nargs="?", type=float, default=0.02)
ap.add_argument("--shipped-json", required=True)
a = ap.parse_args()
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ab", os.path.join(here, "analyze_bench.py"))
ab = importlib.util.module_from_spec(spec)
argv, sys.argv = sys.argv, ["x"]
spec.loader.exec_module(ab)
sys.argv = argv
RANKS = ab.RANKS
POOL = {"species": ["loo"], "genus": ["loo"], "family": ["los"], "order": ["los"], "class": ["los"], "phylum": ["los"]}
ship = json.load(open(a.shipped_json))

d = ab.load([os.path.join(a.bench, f"bench_{m}_s10000.tsv") for m in ("loo", "los", "log")])
d = d[~d["id"].isin(set(open(a.excl).read().split()))].reset_index(drop=True)
ab.NC, ab.NSP, ab.NGE = ab.counts()
d = ab.add_answerable(d)
d["group"] = np.where(d["genus_true"] != "", d["genus_true"], np.where(d["species_true"] != "", d["species_true"], d["id"]))
d["pgroup"] = np.where(d["phylum_true"].isin(["Ascomycota", "Basidiomycota"]), d["phylum_true"], "other")


def frame(dd, rank):
    x = ab.rank_frame(dd[dd["mode"].isin(POOL[rank])], rank)
    return x[x["voted"]]


def fit(x, tgt):
    b = ab.best_threshold(x, tgt)
    return None if b is None else (b["min_ani"], b["min_consistency"])


def apply(x, p):
    return np.zeros(len(x), bool) if p is None else ((x["ani"] >= p[0]) & (x["cons"] >= p[1])).values


def tighter(p, q):
    """The stricter of two (ani, cons) cutoffs, componentwise; an unreachable fit (None) leaves the shipped cutoff."""
    if q is None:
        return p
    return (max(p[0], q[0]), max(p[1], q[1]))


rows = []
groups = d["group"].unique()
for s in range(a.n_splits):
    rng = np.random.default_rng(2000 + s)
    cal_g = set(rng.permutation(groups)[: len(groups) // 2])
    cal, tst = d[d["group"].isin(cal_g)], d[~d["group"].isin(cal_g)]
    for i, rank in enumerate(RANKS):
        xc, xt = frame(cal, rank), frame(tst, rank)
        if len(xc) < 100 or len(xt) < 100:
            continue
        sh = (ship["min_ani"][i], ship["min_consistency"][i])
        mg = fit(xc, a.target + a.margin)
        strat = {"shipped": apply(xt, sh), "fit": apply(xt, fit(xc, a.target)), "margin": apply(xt, mg),
                 "new": apply(xt, tighter(sh, mg))}
        ok = xt["vote_ok"].values
        for name, m in strat.items():
            for grp in ["all", "Ascomycota", "Basidiomycota", "other"]:
                sel = np.ones(len(xt), bool) if grp == "all" else (xt["pgroup"].values == grp)
                rows.append(dict(rank=rank, strategy=name, group=grp, n=int(sel.sum()), called=int((m & sel).sum()),
                                 correct=int(ok[m & sel].sum())))
    print("split", s, "done", file=sys.stderr)

r = pd.DataFrame(rows)
g = r.groupby(["rank", "strategy", "group"]).agg(called=("called", "sum"), correct=("correct", "sum"), n=("n", "sum")).reset_index()
g["precision"] = g["correct"] / g["called"].clip(lower=1)
g["coverage"] = g["called"] / g["n"].clip(lower=1)
g.to_csv(a.out + ".tsv", sep="\t", index=False, float_format="%.4f")

md = [f"# Safety-margin derivation (target {a.target}, margin +{a.margin})\n",
      f"{len(d)} bench queries, {a.n_splits} random genus splits, pool `known`.\n",
      "## Fits on the full `known` pool (in-sample)\n",
      "| rank | shipped (ani, cons) | fit at target | fit at target+margin | new = stricter of shipped and margin fit |", "|---|---|---|---|---|"]
newj = json.loads(json.dumps(ship))
for i, rank in enumerate(RANKS):
    x = frame(d, rank)
    sh = (ship["min_ani"][i], ship["min_consistency"][i])
    f0, f1 = fit(x, a.target), fit(x, a.target + a.margin)
    nw = tighter(sh, f1)
    newj["min_ani"][i], newj["min_consistency"][i] = nw
    fmt = lambda p: "none" if p is None else f"{p[0]:.4f}, {p[1]}"
    md.append(f"| {rank} | {fmt(sh)} | {fmt(f0)} | {fmt(f1)} | {fmt(nw)} |")
json.dump(newj, open(a.out + ".thresholds.json", "w"), indent=1)
for grp in ["all", "Ascomycota", "Basidiomycota", "other"]:
    md.append(f"\n## True phylum group: {grp} (held-out, precision (coverage))\n")
    md.append("| rank | shipped | fit | margin | new |\n|---|---|---|---|---|")
    for rank in RANKS:
        cells = []
        for s in ["shipped", "fit", "margin", "new"]:
            q = g[(g["rank"] == rank) & (g.strategy == s) & (g.group == grp)]
            cells.append("-" if q.empty else f"{q.precision.iloc[0]:.3f} ({q.coverage.iloc[0]:.2f})")
        md.append(f"| {rank} | " + " | ".join(cells) + " |")
open(a.out + ".md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
