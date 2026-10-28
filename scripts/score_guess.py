#!/usr/bin/env python3
"""Score `gmg guess` JSON lines against a truth table.

usage: score_guess.py TRUTH.tsv GUESS.jsonl [GUESS2.jsonl ...] > report.md
TRUTH.tsv: id phylum class order family genus species ... [species_in_bfd genus_in_bfd]
Query id = basename of the query path without extension (.gsk/.fna/.fa/.fastq.gz...).
Reports, per rank: share of queries that got a call (thresholded), precision of those calls,
and the raw top-vote accuracy. Splits by species_in_bfd when that column exists.
"""
import sys, json, os, re
import pandas as pd

RANKS = ["phylum", "class", "order", "family", "genus", "species"]
truth = pd.read_csv(sys.argv[1], sep="\t", dtype=str, keep_default_na=False).set_index("id")
rows = []
for p in sys.argv[2:]:
    for line in open(p):
        d = json.loads(line)
        qid = re.sub(r"(\.fastq\.gz|\.fq\.gz|\.fna|\.fa|\.fasta|\.gsk)+$", "", os.path.basename(d["query"]))
        qid = re.sub(r"_(0\.5|1|2|5|10|20)x.*$", "", qid)  # simulated-read names: <id>_<cov>x
        if qid not in truth.index:
            continue
        t = truth.loc[qid]
        r = {"id": qid, "file": os.path.basename(p), "top_shared": (d["top"] or {}).get("shared", 0), "ani": (d["top"] or {}).get("ani", 0)}
        for k in RANKS:
            c = d["calls"][k]
            r[f"{k}_called"] = c["taxon"] is not None
            r[f"{k}_call_ok"] = c["taxon"] == t[k] if t[k] else None
            r[f"{k}_vote_ok"] = c["vote"] == t[k] if t[k] else None
        for extra in ("species_in_bfd", "genus_in_bfd", "phylum"):
            if extra in t.index:
                r[extra] = t[extra]
        rows.append(r)
df = pd.DataFrame(rows)
print(f"# guess scores: {len(df)} queries from {df['file'].nunique()} file(s)\n")


def table(d):
    out = ["| rank | n | called | precision of calls | top-vote accuracy |", "|---|---|---|---|---|"]
    for k in RANKS:
        x = d[d[f"{k}_vote_ok"].notna()]
        if x.empty:
            continue
        called = x[x[f"{k}_called"]]
        prec = called[f"{k}_call_ok"].astype(float).mean() if len(called) else float("nan")
        out.append(f"| {k} | {len(x)} | {len(called)/len(x):.3f} | {prec:.3f} | {x[f'{k}_vote_ok'].astype(float).mean():.3f} |")
    return "\n".join(out)


for f, g in df.groupby("file"):
    print(f"## {f}\n")
    print(table(g), "\n")
    if "species_in_bfd" in g:
        for flag, gg in g.groupby("species_in_bfd"):
            print(f"species_in_bfd = {flag} (n={len(gg)})\n")
            print(table(gg), "\n")
