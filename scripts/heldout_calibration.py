#!/usr/bin/env python3
"""Held-out calibration (issue #38): do thresholds fitted on some queries keep their precision on others?

usage: heldout_calibration.py BENCH_DIR META.tsv EXCLUDE_IDS.txt OUT_PREFIX [N_SPLITS=20] [TARGET=0.95]
  BENCH_DIR has bench_{loo,los,log}_s10000.tsv from `gmg bench`. Genomes in EXCLUDE_IDS are dropped as queries.
Two held-out designs, for both threshold sets (`known`: reference has relatives; `novel`: query may be a new genus):
  * random splits BY GENUS: queries of half the genera calibrate, the other half test (20 repeats, different seeds);
  * leave-clade-out: calibrate on Ascomycota queries, test on all other phyla.
Reports precision and coverage of the calls the tool would make on the TEST queries, next to the in-sample numbers the
earlier calibration reported, and how stable the fitted thresholds are.
"""
import importlib.util, os, sys, json
import numpy as np
import pandas as pd

bench, meta, excl, out = sys.argv[1:5]
n_splits = int(sys.argv[5]) if len(sys.argv) > 5 else 20
target = float(sys.argv[6]) if len(sys.argv) > 6 else 0.95
here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ab", os.path.join(here, "analyze_bench.py"))
ab = importlib.util.module_from_spec(spec)
sys.argv_backup, sys.argv = sys.argv, ["x"]
spec.loader.exec_module(ab)
sys.argv = sys.argv_backup
RANKS = ab.RANKS
ab.NC, ab.NSP, ab.NGE = ab.counts(meta)

d = ab.load([os.path.join(bench, f"bench_{m}_s10000.tsv") for m in ("loo", "los", "log")])
ex = set(open(excl).read().split())
d = d[~d["id"].isin(ex)].reset_index(drop=True)
d = ab.add_answerable(d)
d["group"] = np.where(d["genus_true"] != "", d["genus_true"], np.where(d["species_true"] != "", d["species_true"], d["id"]))
POOLS = {
    "known": {"species": ["loo"], "genus": ["loo"], "family": ["los"], "order": ["los"], "class": ["los"], "phylum": ["los"]},
    "novel": {"species": ["loo", "los"], "genus": ["los", "log"], "family": ["log"], "order": ["log"], "class": ["log"], "phylum": ["log"]},
}


def pool_frame(dd, rank, modes, seed):
    frames = [ab.rank_frame(dd[dd["mode"] == m], rank) for m in modes]
    frames = [f for f in frames if len(f)]
    if not frames:
        return pd.DataFrame()
    if len(frames) > 1:
        n = min(len(f) for f in frames)
        frames = [f.sample(n=n, random_state=seed) for f in frames]
    return pd.concat(frames)


def evaluate(x, ani, cons):
    v = x[x["voted"]]
    if v.empty:
        return np.nan, np.nan, 0
    sel = v[(v["cons"] >= cons) & (v["ani"] >= ani)]
    return (sel["vote_ok"].mean() if len(sel) else np.nan), len(sel) / len(v), len(sel)


def fit(dd, rank, pool, seed):
    return ab.best_threshold(pool_frame(dd, rank, POOLS[pool][rank], seed), target)


rows, thr_rows = [], []
groups = d["group"].unique()
for pool in POOLS:
    for rank in RANKS:
        # in-sample reference: fit and evaluate on everything
        full = pool_frame(d, rank, POOLS[pool][rank], 1)
        b_all = ab.best_threshold(full, target) if len(full) else None
        ins = evaluate(full, b_all["min_ani"], b_all["min_consistency"]) if b_all else (np.nan, np.nan, 0)
        res = []
        for s in range(n_splits):
            rng = np.random.default_rng(1000 + s)
            cal_groups = set(rng.permutation(groups)[: len(groups) // 2])
            cal = d[d["group"].isin(cal_groups)]
            tst = d[~d["group"].isin(cal_groups)]
            b = fit(cal, rank, pool, s)
            if not b:
                res.append((np.nan, np.nan, 0, np.nan, np.nan)); continue
            p, c, n = evaluate(pool_frame(tst, rank, POOLS[pool][rank], s), b["min_ani"], b["min_consistency"])
            res.append((p, c, n, b["min_ani"], b["min_consistency"]))
        r = np.array(res, dtype=float)
        # leave-clade-out
        asc = d[d["phylum_true"] == "Ascomycota"]; oth = d[d["phylum_true"] != "Ascomycota"]
        bc = fit(asc, rank, pool, 1)
        pc = evaluate(pool_frame(oth, rank, POOLS[pool][rank], 1), bc["min_ani"], bc["min_consistency"]) if bc else (np.nan, np.nan, 0)
        rows.append(dict(pool=pool, rank=rank, insample_prec=ins[0], insample_cov=ins[1],
                         heldout_prec=np.nanmean(r[:, 0]), prec_lo=np.nanpercentile(r[:, 0], 2.5), prec_hi=np.nanpercentile(r[:, 0], 97.5),
                         heldout_cov=np.nanmean(r[:, 1]), cov_lo=np.nanpercentile(r[:, 1], 2.5), cov_hi=np.nanpercentile(r[:, 1], 97.5),
                         min_n=int(np.nanmin(r[:, 2])) if np.isfinite(r[:, 2]).any() else 0, fits=int(np.isfinite(r[:, 3]).sum()),
                         ani_mean=np.nanmean(r[:, 3]), ani_sd=np.nanstd(r[:, 3]),
                         clade_prec=pc[0], clade_cov=pc[1], clade_n=pc[2]))
df = pd.DataFrame(rows)
df.to_csv(out + ".tsv", sep="\t", index=False, float_format="%.4f")

fmt = lambda x: "-" if pd.isna(x) else f"{x:.3f}"
md = [f"# Held-out calibration (target precision {target})\n",
      f"{len(d)} bench queries (clean labels: {len(ex)} suspect genomes dropped), {len(groups)} genus groups, {n_splits} random splits by genus (half the genera calibrate, the other half test).\n",
      "`held-out precision` = precision of the calls made with thresholds fitted on the other half; brackets = 2.5-97.5 percentile over the splits. `clade-out` = fit on Ascomycota queries, test on all other phyla.\n"]
for pool in POOLS:
    md.append(f"\n## `{pool}` thresholds\n")
    md.append("| rank | in-sample precision | held-out precision | held-out coverage | fitted min ANI (mean ± sd) | clade-out precision (coverage, n) |\n|---|---|---|---|---|---|")
    for _, r in df[df.pool == pool].iterrows():
        md.append(f"| {r['rank']} | {fmt(r.insample_prec)} | {fmt(r.heldout_prec)} [{fmt(r.prec_lo)}, {fmt(r.prec_hi)}] | {fmt(r.heldout_cov)} [{fmt(r.cov_lo)}, {fmt(r.cov_hi)}] | "
                  f"{fmt(r.ani_mean)} ± {fmt(r.ani_sd)} | {fmt(r.clade_prec)} ({fmt(r.clade_cov)}, {int(r.clade_n)}) |")
md.append("\nReading: a held-out precision at or above the target means the thresholds generalise to genera they were not fitted on; below it means the in-sample numbers were optimistic. "
          "Ranks with no fitted threshold at the target (`-`) never reach it for that pool.\n")
open(out + ".md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
