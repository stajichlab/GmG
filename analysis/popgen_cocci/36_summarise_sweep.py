#!/usr/bin/env python3
"""36_summarise_sweep.py : tabulate the threshold sweep (35_threshold_sweep.sbatch): per setting, how many samples are flagged
or excluded, the cut verdict, and how the flagged set differs from the default.   usage: 36_summarise_sweep.py cocci|bd"""
import sys, json, glob, os
import pandas as pd
R = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsketch_runs/"
ds = sys.argv[1]
base = None
rows = []
names = sorted(os.listdir(R + f"sweep_{ds}"), key=lambda n: (n != "base", n.rsplit("_", 1)[0], float(n.rsplit("_", 1)[1]) if "_" in n else 0))
for n in names:
    d = R + f"sweep_{ds}/{n}/"
    if not os.path.exists(d + "qc.tsv"):
        rows.append(dict(setting=n, note="no output")); continue
    q = pd.read_csv(d + "qc.tsv", sep="\t").fillna("")
    fl = {(i, f) for i, s in zip(q["id"], q["flags"]) for f in str(s).split(",") if f}
    ex = set(q["id"][q["status"] == "excluded"])
    j = json.load(open(d + "popsketch.json"))
    g = (j.get("clusters") or {}).get("gap") or {}
    if n == "base": base = (fl, ex)
    rows.append(dict(setting=n, flagged_samples=len({i for i, _ in fl}), excluded=len(ex), verdict=g.get("verdict", "-"),
                     gap_over_floor=round(g.get("gap_over_floor", 0), 1) if g else "-",
                     flags_added=len(fl - base[0]) if base else "", flags_lost=len(base[0] - fl) if base else "",
                     excl_added=len(ex - base[1]) if base else "", excl_lost=len(base[1] - ex) if base else ""))
print(pd.DataFrame(rows).to_string(index=False))
