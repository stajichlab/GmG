#!/usr/bin/env python3
"""Score stress-test results (issue #41): usage stress_score.py RESULTS.tsv TRUTH.tsv OUT.md"""
import sys
import pandas as pd

res = pd.read_csv(sys.argv[1], sep="\t", dtype=str, keep_default_na=False)
tr = pd.read_csv(sys.argv[2], sep="\t", dtype=str, keep_default_na=False)
d = tr.merge(res, left_on="id", right_on="query_id")
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
md = ["# Stress tests: degraded known genomes\n",
      f"{len(d)} inputs from {tr[tr.host != ''].host.nunique()} chromosome-level assemblies (genus present in BFD, not in the index). Truth = the host lineage. "
      "Each cell: correct/called/total for genus, then phylum.\n"]
md.append("| variant | param | n | phylum called (correct) | genus called (correct) | species called (correct) | median hashes | flags seen |\n|---|---|---|---|---|---|---|---|")
for (var, par), g in d[d.variant != "negative"].groupby(["variant", "param"], sort=False):
    def cell(r):
        c = g[g[f"call_{r}"] != ""]; ok = (c[f"call_{r}"] == c[r]).sum()
        return f"{len(c)}/{len(g)} ({ok})"
    flags = g["qc_flags"].str.split(",").explode().replace("", pd.NA).value_counts().head(3)
    md.append(f"| {var} | {par} | {len(g)} | {cell('phylum')} | {cell('genus')} | {cell('species')} | {int(pd.to_numeric(g.n_hashes).median())} | {', '.join(f'{k}×{v}' for k, v in flags.items())} |")
neg = d[d.variant == "negative"]
md.append("\n## Negative controls\n")
for _, r in neg.iterrows():
    md.append(f"* `{r.id}`: status `{r.status}`, phylum call `{r.call_phylum or '-'}`, flags `{r.qc_flags}`")
md.append("\nA tentative genus is not a call; `call_*` columns only. Wrong calls are those where the called taxon differs from the host lineage; for the contamination series the host is the majority genome.\n")
open(sys.argv[3], "w").write("\n".join(md) + "\n")
print("\n".join(md))
