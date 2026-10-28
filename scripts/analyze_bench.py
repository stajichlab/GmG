#!/usr/bin/env python3
"""Summarise and calibrate `gmg bench` output.

usage: analyze_bench.py OUT_PREFIX BENCH.tsv [BENCH2.tsv ...]
  Each TSV may be one mode (loo/los/log); the mode is taken from the file name (bench_<mode>*.tsv).
Writes OUT_PREFIX.md (tables), OUT_PREFIX.thresholds.json (calibrated for `gmg guess --thresholds`).

For every rank we ask: when the vote names a taxon and the evidence is >= (shared, consistency),
how often is the vote right? Truth is known only when the rank is filled in the metadata.
Thresholds target PRECISION (default 0.95) while keeping as much coverage as possible.
"""
import sys, os, re, json
import pandas as pd
import numpy as np

RANKS = ["phylum", "class", "order", "family", "genus", "species"]
TARGET = float(os.environ.get("GMG_PRECISION", "0.95"))


def load(paths):
    frames = []
    for p in paths:
        m = re.search(r"bench_(loo|los|log)", os.path.basename(p))
        d = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        for c in ["query_hashes", "top_shared", "n_ref"]:
            d[c] = pd.to_numeric(d[c])
        for c in ["containment_query", "containment_ref", "ani"] + [f"{r}_cons" for r in RANKS]:
            d[c] = pd.to_numeric(d[c])
        d["mode"] = m.group(1) if m else "loo"
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def rank_frame(d, r):
    x = d[d[f"{r}_true"] != ""].copy()
    x["vote_ok"] = (x[f"{r}_vote"] == x[f"{r}_true"])
    x["voted"] = x[f"{r}_vote"] != ""
    x["cons"] = x[f"{r}_cons"]
    return x


def best_threshold(x, target):
    """Pick (min_shared, min_consistency) maximising coverage subject to precision >= target."""
    x = x[x["voted"]]
    if x.empty:
        return None
    best = None
    shared_grid = sorted(set(np.unique(np.percentile(x["top_shared"], np.arange(0, 100, 2))).astype(int)) | {1, 5, 10, 20, 50, 100})
    for cons in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        for sh in shared_grid:
            sel = x[(x["top_shared"] >= sh) & (x["cons"] >= cons)]
            if len(sel) < 20:
                continue
            prec = sel["vote_ok"].mean()
            cov = len(sel) / len(x)
            if prec >= target and (best is None or cov > best["coverage"]):
                best = dict(min_shared=int(sh), min_consistency=cons, precision=float(prec), coverage=float(cov), n=int(len(sel)))
    return best


def main():
    out = sys.argv[1]
    d = load(sys.argv[2:])
    md = [f"# GmG benchmark summary\n\nqueries: {len(d)}; modes: {sorted(d['mode'].unique())}; precision target {TARGET}\n"]
    thr = {}
    for mode in sorted(d["mode"].unique()):
        dm = d[d["mode"] == mode]
        md.append(f"\n## mode `{mode}` (n={len(dm)})\n")
        md.append("| rank | n with truth | top-hit correct | vote correct | no hit |\n|---|---|---|---|---|")
        for r in RANKS:
            x = rank_frame(dm, r)
            if x.empty:
                continue
            top_ok = (x[f"{r}_top"] == x[f"{r}_true"]).mean()
            md.append(f"| {r} | {len(x)} | {top_ok:.3f} | {x['vote_ok'].mean():.3f} | {(x['top_shared']==0).mean():.3f} |")
        md.append("\nby phylum (vote correct):\n")
        md.append("| phylum | n | " + " | ".join(RANKS) + " |\n|---|---|" + "---|" * 6)
        for ph, g in dm.groupby("phylum_true"):
            if len(g) < 20 or ph == "":
                continue
            cells = []
            for r in RANKS:
                x = rank_frame(g, r)
                cells.append(f"{x['vote_ok'].mean():.3f}" if len(x) >= 10 else "-")
            md.append(f"| {ph} | {len(g)} | " + " | ".join(cells) + " |")
    # calibration from the matching mode per rank: species <- los is too harsh; use loo for species/genus,
    # los for genus when available, log for order/family/class/phylum when available
    pref = {"species": ["loo"], "genus": ["los", "loo"], "family": ["log", "los", "loo"],
            "order": ["log", "los", "loo"], "class": ["log", "los", "loo"], "phylum": ["log", "los", "loo"]}
    th = dict(min_shared=[5] * 6, min_ani=[0.0] * 6, min_consistency=[0.6] * 6, margin=0.9, top_n=10)
    md.append(f"\n## calibration (precision >= {TARGET})\n\n| rank | basis mode | min_shared | min_consistency | precision | coverage | n |\n|---|---|---|---|---|---|---|")
    for i, r in enumerate(RANKS):
        for mode in pref[r]:
            dm = d[d["mode"] == mode]
            if dm.empty:
                continue
            b = best_threshold(rank_frame(dm, r), TARGET)
            if b:
                th["min_shared"][i] = b["min_shared"]
                th["min_consistency"][i] = b["min_consistency"]
                md.append(f"| {r} | {mode} | {b['min_shared']} | {b['min_consistency']} | {b['precision']:.3f} | {b['coverage']:.3f} | {b['n']} |")
                break
        else:
            md.append(f"| {r} | - | (default) | (default) | - | - | - |")
    open(out + ".md", "w").write("\n".join(md) + "\n")
    json.dump(th, open(out + ".thresholds.json", "w"), indent=1)
    print("\n".join(md))


main()
