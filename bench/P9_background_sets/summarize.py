#!/usr/bin/env python3
"""P9 summary from out/work/{profile,results,exits}.tsv: outcome per run (latest attempt), per-set counts, and median/max of each fetch stage for the ok runs.
usage: summarize.py WORKDIR MANIFEST [SETCSV ...]   (a SETCSV has a run_accession column; the set name is the file name without .csv)"""
import csv, statistics, sys, os
work, manifest, sets = sys.argv[1], sys.argv[2], sys.argv[3:]
prof = {}
for r in csv.DictReader(open(f"{work}/profile.tsv"), delimiter="\t"): prof[r["run"]] = r  # later attempts replace earlier ones
exits = {}
for l in open(f"{work}/exits.tsv"):
    p = l.rstrip("\n").split("\t")
    if len(p) >= 2: exits[p[0]] = p[1]
runs = [r.get("run_accession") or r.get("run") for r in csv.DictReader(open(manifest))]
print(f"manifest {len(runs)} runs; latest outcome: " + ", ".join(f"{k} {v}" for k, v in sorted(__import__('collections').Counter(exits.get(r, 'none') for r in runs).items())))
for s in sets:
    ids = {r["run_accession"] for r in csv.DictReader(open(s))}
    print(f"set {os.path.basename(s)[:-4]}: {len(ids)} runs, ok {sum(exits.get(i) == '0' for i in ids)}")
ok = [prof[r] for r in runs if exits.get(r) == "0" and r in prof]
print(f"stage seconds over {len(ok)} ok runs (median, max):")
for c in ("seconds", "transfer_s", "convert_s", "count_s", "cut_s", "hook_s", "tree_peak_mib"):
    v = [float(r[c]) for r in ok if r.get(c)]
    if v: print(f"  {c}: {statistics.median(v):.1f}, {max(v):.1f}")
print("sampling: " + ", ".join(f"{k} {v}" for k, v in sorted(__import__('collections').Counter(r["sampling"] for r in ok).items())))
