#!/usr/bin/env python3
"""Summarise and calibrate `gmg bench` output.

usage: analyze_bench.py OUT_PREFIX BENCH.tsv [BENCH2.tsv ...]     (env GMG_META = meta.tsv of the reference set)
  The mode (loo/los/log) is taken from the file name (bench_<mode>*.tsv).
Writes OUT_PREFIX.md and OUT_PREFIX.thresholds.json (for `gmg guess --thresholds`).

"Answerable": in leave-one-out a rank is answerable only if another genome carries the true taxon;
in leave-species-out only if a genome of a *different* species does; in leave-genus-out only if a genome
of a different genus does. Unanswerable queries are counted as errors in `all` columns and left out of
the `answerable` columns, so a singleton species cannot make the tool look worse than it is.

Calibration targets PRECISION (GMG_PRECISION, default 0.95) of the calls the tool would make, using
features that do not depend on `scaled`: ANI (containment ANI of the top hit) and vote consistency.
Set GMG_EXCLUDE=file (ids, one per line) to drop suspected mislabeled genomes from the queries.
Pools: species = loo + los (a novel species must not get a species call); genus = los + log;
family/order/class/phylum = log.
"""
import sys, os, re, json
import pandas as pd
import numpy as np

RANKS = ["phylum", "class", "order", "family", "genus", "species"]
TARGET = float(os.environ.get("GMG_PRECISION", "0.95"))
META = os.environ.get("GMG_META", "/bigdata/stajichlab/shared/projects/GmG_runs/index/meta.tsv")


def load(paths):
    frames = []
    for p in paths:
        m = re.search(r"bench_(loo|los|log)", os.path.basename(p))
        d = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        for c in ["query_hashes", "top_shared", "n_ref"]:
            d[c] = pd.to_numeric(d[c])
        for c in ["containment_query", "containment_ref", "ani"] + [f"{r}_cons" for r in RANKS]:
            d[c] = pd.to_numeric(d[c])
        d["ani"] = d["containment_query"] ** (1.0 / 21)   # query-containment ANI (k=21); the TSV column may use the older max() form
        d["mode"] = m.group(1) if m else "loo"
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def counts():
    meta = pd.read_csv(META, sep="\t", dtype=str, keep_default_na=False)
    cols = dict(zip(RANKS, ["phylum", "class", "order", "family", "genus", "species"]))
    n = {r: meta[cols[r]].value_counts().to_dict() for r in RANKS}
    n_sp = {r: meta.groupby([cols[r], "species"]).size().to_dict() for r in RANKS}
    n_ge = {r: meta.groupby([cols[r], "genus"]).size().to_dict() for r in RANKS}
    return n, n_sp, n_ge


NC, NSP, NGE = counts()
EXCLUDE = set(open(os.environ['GMG_EXCLUDE']).read().split()) if os.environ.get('GMG_EXCLUDE') else set()


def add_answerable(d):
    mode, sp, ge = d["mode"].values, d["species_true"].values, d["genus_true"].values
    for r in RANKS:
        tt = d[f"{r}_true"].values
        out = np.zeros(len(d), dtype=bool)
        for i, t in enumerate(tt):
            if not t:
                continue
            m = mode[i]
            if m == "loo":
                out[i] = NC[r].get(t, 0) >= 2
            elif m == "los":
                out[i] = r != "species" and NC[r].get(t, 0) - NSP[r].get((t, sp[i]), 0) >= 1
            else:
                out[i] = r not in ("genus", "species") and NC[r].get(t, 0) - NGE[r].get((t, ge[i]), 0) >= 1
        d[f"{r}_ans"] = out
    return d


def rank_frame(d, r):
    x = d[d[f"{r}_true"] != ""].copy()
    x["vote_ok"] = (x[f"{r}_vote"] == x[f"{r}_true"])
    x["top_ok"] = (x[f"{r}_top"] == x[f"{r}_true"])
    x["voted"] = x[f"{r}_vote"] != ""
    x["cons"] = x[f"{r}_cons"]
    x["answerable"] = x[f"{r}_ans"]
    return x


def best_threshold(x, target):
    """(min_ani, min_consistency) maximising coverage subject to precision >= target."""
    x = x[x["voted"]]
    if len(x) < 50:
        return None
    best = None
    ani_grid = np.unique(np.round(np.percentile(x["ani"], np.arange(0, 100, 1.0)), 4))
    for cons in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        xc = x[x["cons"] >= cons]
        for a in ani_grid:
            sel = xc[xc["ani"] >= a]
            if len(sel) < 30:
                break
            prec = sel["vote_ok"].mean()
            cov = len(sel) / len(x)
            if prec >= target and (best is None or cov > best["coverage"]):
                best = dict(min_ani=float(a), min_consistency=cons, precision=float(prec), coverage=float(cov), n=int(len(sel)))
    return best


def main():
    out = sys.argv[1]
    d = load(sys.argv[2:])
    if EXCLUDE:
        d = d[~d['id'].isin(EXCLUDE)].reset_index(drop=True)
    d = add_answerable(d)
    md = [f"# GmG benchmark summary\n\nqueries: {len(d)}; modes: {sorted(d['mode'].unique())}; precision target {TARGET}\n"]
    for mode in sorted(d["mode"].unique()):
        dm = d[d["mode"] == mode]
        md.append(f"\n## mode `{mode}` (n={len(dm)})\n")
        md.append("| rank | n | answerable | vote correct (all) | vote correct (answerable) | top-hit correct (answerable) |\n|---|---|---|---|---|---|")
        for r in RANKS:
            x = rank_frame(dm, r)
            if x.empty:
                continue
            a = x[x["answerable"]]
            if a.empty:
                md.append(f"| {r} | {len(x)} | 0 | {x['vote_ok'].mean():.3f} | - | - |")
            else:
                md.append(f"| {r} | {len(x)} | {len(a)} ({len(a)/len(x):.2f}) | {x['vote_ok'].mean():.3f} | {a['vote_ok'].mean():.3f} | {a['top_ok'].mean():.3f} |")
        md.append("\nby phylum (vote correct, answerable only; n shown for the rank-species column):\n")
        md.append("| phylum | n | " + " | ".join(RANKS) + " |\n|---|---|" + "---|" * 6)
        for ph, g in dm.groupby("phylum_true"):
            if len(g) < 20 or ph == "":
                continue
            cells = []
            for r in RANKS:
                x = rank_frame(g, r)
                a = x[x["answerable"]]
                cells.append(f"{a['vote_ok'].mean():.3f} (n={len(a)})" if len(a) >= 10 else "-")
            md.append(f"| {ph} | {len(g)} | " + " | ".join(cells) + " |")
    # GMG_POOLS=known : the reference holds close relatives (species/genus from leave-one-out, higher ranks from leave-species-out)
    # GMG_POOLS=novel : the query may be a new species/genus (species from loo+los, genus from los+log, higher ranks from leave-genus-out)
    if os.environ.get("GMG_POOLS", "novel") == "known":
        pools = {"species": ["loo"], "genus": ["loo"], "family": ["los"], "order": ["los"], "class": ["los"], "phylum": ["los"]}
    else:
        pools = {"species": ["loo", "los"], "genus": ["los", "log"], "family": ["log"], "order": ["log"], "class": ["log"], "phylum": ["log"]}
    th = dict(min_shared=[5] * 6, min_ani=[0.0] * 6, min_consistency=[0.6] * 6, margin=0.9, top_n=10)
    md.append(f"\n## calibration (precision >= {TARGET}, counting unanswerable queries as errors)\n\n| rank | pool | min_ani | min_consistency | precision | coverage | n |\n|---|---|---|---|---|---|---|")
    for i, r in enumerate(RANKS):
        frames = [rank_frame(d[d["mode"] == m], r) for m in pools[r] if (d["mode"] == m).any()]
        if not frames:
            continue
        # equal weight per mode when pooling
        if len(frames) > 1:
            n = min(len(f) for f in frames)
            frames = [f.sample(n=n, random_state=1) for f in frames]
        x = pd.concat(frames)
        b = best_threshold(x, TARGET)
        if b:
            th["min_ani"][i] = b["min_ani"]
            th["min_consistency"][i] = b["min_consistency"]
            md.append(f"| {r} | {'+'.join(pools[r])} | {b['min_ani']:.4f} | {b['min_consistency']} | {b['precision']:.3f} | {b['coverage']:.3f} | {b['n']} |")
        else:
            md.append(f"| {r} | {'+'.join(pools[r])} | (no threshold reaches the target) | | | | |")
    open(out + ".md", "w").write("\n".join(md) + "\n")
    json.dump(th, open(out + ".thresholds.json", "w"), indent=1)
    print("\n".join(md))


main()
