#!/usr/bin/env python3
"""Score a `gmg run` results.tsv against a truth table (issue #25 groundwork; used by every validation dataset).

usage: validate_run.py RESULTS.tsv TRUTH.tsv OUT_PREFIX [--by COL ...] [--title TEXT]
TRUTH.tsv: id + any of phylum class order family genus species (empty = unknown, not scored) + optional grouping columns.
For every rank: n with truth, share of inputs given a call, PRECISION of the calls made, recall (correct calls / n with truth),
and the raw top-vote accuracy. Breakdowns: for each --by column (default: species_in_bfd genus_in_bfd level phylum) and for
quality bins (total_len, n50) when present in the truth table. Also: status counts, and for known/novel taxa the false-call rate
("called a species/genus that is not the truth").
Also: confusion lists (most frequent wrong genus/species calls), calibration (precision by ANI and vote consistency) and, for reads, accuracy by amount of sequence.
Writes OUT_PREFIX.md and OUT_PREFIX.tsv (long format: group, rank, n, called, precision, recall, vote_acc).
"""
import sys, argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("results"); ap.add_argument("truth"); ap.add_argument("out")
ap.add_argument("--by", nargs="*", default=["species_in_bfd", "genus_in_bfd", "level", "phylum"])
ap.add_argument("--title", default="Validation")
a = ap.parse_args()
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
res = pd.read_csv(a.results, sep="\t", dtype=str, keep_default_na=False)
tru = pd.read_csv(a.truth, sep="\t", dtype=str, keep_default_na=False)
d = res.merge(tru, left_on="query_id", right_on="id", how="inner", suffixes=("", "_t"))
unmatched = len(res) - len(d)
for r in RANKS:
    d[f"ok_{r}"] = d[f"call_{r}"] == d[r]
    d[f"has_{r}"] = d[r] != ""
    d[f"called_{r}"] = d[f"call_{r}"] != ""
for c in ("total_len", "n50"):
    if c in d: d[c] = pd.to_numeric(d[c], errors="coerce")
if "n50" in d:
    d["n50_bin"] = pd.cut(d["n50"], [-1, 10_000, 100_000, 1_000_000, 1e12], labels=["<10 kb", "10-100 kb", "0.1-1 Mb", ">= 1 Mb"]).astype(str)
if "total_len" in d:
    d["size_bin"] = pd.cut(d["total_len"], [-1, 10e6, 25e6, 50e6, 1e12], labels=["<10 Mb", "10-25 Mb", "25-50 Mb", ">= 50 Mb"]).astype(str)
d["conf"] = d["status"]
scored = d[d["status"] != "error"]


def metrics(g, rank):
    t = g[g[f"has_{rank}"]]
    if t.empty: return None
    called = t[t[f"called_{rank}"]]
    prec = called[f"ok_{rank}"].mean() if len(called) else np.nan
    return dict(n=len(t), called=len(called) / len(t), precision=prec, recall=t[f"ok_{rank}"].mean())


def table(g, label):
    rows = []
    for r in RANKS:
        m = metrics(g, r)
        if m: rows.append(dict(group=label, rank=r, **m))
    return rows


def wilson(p, n, z=1.96):
    if n == 0 or pd.isna(p): return (np.nan, np.nan)
    c = (p + z * z / (2 * n)) / (1 + z * z / n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (c - h, c + h)


def md_table(rows):
    out = ["| rank | n | called | precision [95% CI] | recall |", "|---|---|---|---|---|"]
    for r in rows:
        n_called = round(r["called"] * r["n"])
        lo, hi = wilson(r["precision"], n_called)
        p = "-" if pd.isna(r["precision"]) else f"{r['precision']:.3f} [{lo:.3f}, {hi:.3f}]"
        out.append(f"| {r['rank']} | {r['n']} | {r['called']:.3f} | {p} | {r['recall']:.3f} |")
    return "\n".join(out)


allrows = table(scored, "all")
md = [f"# {a.title}\n", f"{len(res)} inputs scored ({unmatched} without a truth row); status: {dict(res['status'].value_counts())}\n",
      "Precision = correct among the calls made; called = share of inputs with a call at that rank; recall = correct calls / inputs with truth. "
      "Errors are excluded. Truth = the lineage in the truth table (NCBI taxonomy or curated); where it is blank the rank is not scored.\n", "## All inputs\n", md_table(allrows)]
for col in a.by + [c for c in ("n50_bin", "size_bin") if c in d]:
    if col not in d: continue
    md.append(f"\n## By `{col}`\n")
    for val, g in scored.groupby(col):
        rows = table(g, f"{col}={val}")
        if rows and len(g) >= 5:
            md.append(f"\n**{col} = {val}** (n = {len(g)})\n\n" + md_table(rows))
        allrows += rows
# false-call rate for novel taxa
if "species_in_bfd" in d:
    nov = scored[(scored["species_in_bfd"] == "False") & (scored["species"] != "")]
    if len(nov):
        md.append(f"\n## Named species NOT in the reference (n = {len(nov)})\n\nSpecies called for {int(nov['called_species'].sum())} of {len(nov)} (every one of those is wrong by construction: the species is not in the reference); "
                  f"genus called for {int(nov['called_genus'].sum())}, correct {int((nov['called_genus'] & nov['ok_genus']).sum())}.\n")
if "genus_in_bfd" in d:
    ng = scored[(scored["genus_in_bfd"] == "False") & (scored["genus"] != "")]
    if len(ng):
        md.append(f"\n## Genus NOT in the reference (n = {len(ng)})\n\nGenus called for {int(ng['called_genus'].sum())} (wrong by construction); family called {int(ng['called_family'].sum())}, correct {int((ng['called_family'] & ng['ok_family']).sum())}; "
                  f"order called {int(ng['called_order'].sum())}, correct {int((ng['called_order'] & ng['ok_order']).sum())}; phylum called {int(ng['called_phylum'].sum())}, correct {int((ng['called_phylum'] & ng['ok_phylum']).sum())}.\n")
# confusion: the most frequent wrong calls at genus and species (truth -> call)
for r in ("genus", "species"):
    w = scored[scored[f"has_{r}"] & scored[f"called_{r}"] & ~scored[f"ok_{r}"]]
    if len(w):
        c = w.groupby([r, f"call_{r}"]).size().sort_values(ascending=False).head(15)
        md.append(f"\n## Most frequent wrong {r} calls ({len(w)} wrong of {int((scored[f'has_{r}'] & scored[f'called_{r}']).sum())} called)\n\n| truth | called | n |\n|---|---|---|")
        md += [f"| {t} | {k} | {n} |" for (t, k), n in c.items()]
# calibration: precision of genus/species calls by the reported ANI and by the vote consistency
for r in ("genus", "species"):
    t = scored[scored[f"has_{r}"] & scored[f"called_{r}"]].copy()
    for col, bins, lab in (("ani", [0, .90, .95, .98, .995, 1.01], "top-hit ANI"), (f"conf_{r}", [-1, .6, .8, .9, .95, 1.01], f"{r} vote consistency")):
        if col in t and len(t):
            t["_v"] = pd.to_numeric(t[col], errors="coerce")
            t["_b"] = pd.cut(t["_v"], bins, right=False).astype(str)
            g = t.groupby("_b").agg(n=(f"ok_{r}", "size"), precision=(f"ok_{r}", "mean"))
            md.append(f"\n## Calibration: {r} precision by {lab}\n\n| {lab} | n | precision |\n|---|---|---|")
            md += [f"| {b} | {int(x.n)} | {x.precision:.3f} |" for b, x in g.iterrows()]
            allrows += [dict(group=f"calib_{col}={b}", rank=r, n=int(x.n), called=1.0, precision=x.precision, recall=np.nan) for b, x in g.iterrows()]
# reads: accuracy against sequencing depth
if "kind" in d and (scored["kind"] == "reads").any() and "total_bases" in d:
    rd = scored[scored["kind"] == "reads"].copy()
    rd["_tb"] = pd.to_numeric(rd["total_bases"], errors="coerce")
    rd["_bin"] = pd.cut(rd["_tb"], [-1, 1e7, 1e8, 1e9, 1e12], labels=["<10 Mb", "10-100 Mb", "0.1-1 Gb", ">= 1 Gb"]).astype(str)
    for val, g in rd.groupby("_bin"):
        rows = table(g, f"reads_bases={val}")
        if rows and len(g) >= 5:
            md.append(f"\n**reads, {val} of sequence** (n = {len(g)})\n\n" + md_table(rows))
        allrows += rows
pd.DataFrame(allrows).to_csv(a.out + ".tsv", sep="\t", index=False, float_format="%.4f")
open(a.out + ".md", "w").write("\n".join(md) + "\n")
print("\n".join(md))
