#!/usr/bin/env python3
"""Score read-based guesses by coverage / depth.  usage: score_reads.py TRUTH.tsv KIND GUESS.jsonl [...]
KIND = sim  (query names <id>_<cov>x.gsk)  or  cocci (query dir name carries the read count, e.g. cocci_500000)
Prints, per input file and per coverage (or read count), the fraction of queries given a call and the precision of the
calls, plus raw top-vote accuracy, for phylum/order/genus/species."""
import sys, json, os, re
import pandas as pd
truth = pd.read_csv(sys.argv[1], sep="\t", dtype=str, keep_default_na=False).set_index("id")
kind = sys.argv[2]
RK = ["phylum", "order", "genus", "species"]
rows = []
for p in sys.argv[3:]:
    for line in open(p):
        d = json.loads(line)
        base = os.path.basename(d["query"]); stem = re.sub(r"\.gsk$", "", base)
        if kind == "sim":
            m = re.match(r"(.+)_(0\.5|1|5|20)x$", stem); qid, grp = m.group(1), m.group(2) + "x"
        else:
            qid = stem; grp = re.search(r"cocci_(\d+)", d["query"]).group(1) + " reads"
        if qid not in truth.index: continue
        t = truth.loc[qid]
        r = {"file": os.path.basename(p), "grp": grp, "qh": d["query_hashes"], "shared": (d["top"] or {}).get("shared", 0)}
        for k in RK:
            c = d["calls"][k]
            r[k + "_called"] = c["taxon"] is not None
            r[k + "_ok"] = (c["taxon"] == t[k]) if c["taxon"] is not None and t[k] else None
            r[k + "_vote"] = (c["vote"] == t[k]) if t[k] else None
        if "species_in_bfd" in t.index: r["known"] = t["species_in_bfd"]
        rows.append(r)
df = pd.DataFrame(rows)
for f, g in df.groupby("file"):
    print(f"\n## {f}\n")
    print("| group | n | " + " | ".join(f"{k}: called / precision / top-vote" for k in RK) + " |\n|---|---|" + "---|" * len(RK))
    order = sorted(g["grp"].unique(), key=lambda s: float(re.match(r"[\d.]+", s).group()))
    for grp in order:
        x = g[g["grp"] == grp]
        cells = []
        for k in RK:
            vx = x[x[k + "_vote"].notna()]
            called = x[x[k + "_called"] & x[k + "_ok"].notna()]
            prec = called[k + "_ok"].astype(float).mean() if len(called) else float("nan")
            cells.append(f"{len(called)/max(len(vx),1):.2f} / {prec:.2f} / {vx[k+'_vote'].astype(float).mean():.2f}" if len(vx) else "-")
        print(f"| {grp} | {len(x)} | " + " | ".join(cells) + " |")
