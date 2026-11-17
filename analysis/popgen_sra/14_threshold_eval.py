#!/usr/bin/env python3
"""#159 step 3: evaluate the peak-window share sweep against the assembly size (ground truth for size).
usage (from popsra): 14_threshold_eval.py SLUG [SLUG ...] ; reads thr_<t>/<slug>/qc.tsv and <slug>/runs.tsv
A call is right when est/assembly is in 0.75-1.33. Reports, per threshold: peaks (basis == peak), right among all runs that have a size, and the number of runs whose basis differs from the 0.30 run.
Held-out check: the threshold is chosen on the species in TUNE (default the first half of the list) and scored on the rest; the pick is the lowest threshold within one run of the best right-count (ties go to the one nearest 0.30)."""
import csv, os, sys
slugs = sys.argv[1:]; thr = ["0.20", "0.25", "0.30", "0.35", "0.40"]
truth = {s: {r["run_accession"]: float(r["genome_mb"]) for r in csv.DictReader(open(f"{s}/runs.tsv"), delimiter="\t") if r.get("genome_mb")} for s in slugs}
def load(t, s): return {q["id"]: q for q in csv.DictReader(open(f"thr_{t}/{s}/qc.tsv"), delimiter="\t")}
data = {(t, s): load(t, s) for t in thr for s in slugs}
def score(t, group):
    right = n = peaks = 0
    for s in group:
        for i, q in data[(t, s)].items():
            if i in truth[s] and q["est_size_mb"]:
                n += 1; peaks += q["coverage_basis"] == "peak"; right += 0.75 <= float(q["est_size_mb"]) / truth[s][i] <= 1.33
    return right, n, peaks
def changed(t, group):
    return sum(1 for s in group for i, q in data[(t, s)].items() if i in data[("0.30", s)] and q["coverage_basis"] != data[("0.30", s)][i]["coverage_basis"])
half = len(slugs) // 2
for name, group in (("all", slugs), ("tune", slugs[:half]), ("held-out", slugs[half:])):
    print(f"== {name}: {', '.join(group)}")
    for t in thr:
        r, n, p = score(t, group); print(f"  {t}: right {r}/{n} ({100*r/n:.1f}%), peaks {p}, basis differs from 0.30: {changed(t, group)}")
best = max(score(t, slugs[:half])[0] for t in thr)
pick = min((t for t in thr if score(t, slugs[:half])[0] >= best - 1), key=lambda t: abs(float(t) - 0.3))
r, n, _ = score(pick, slugs[half:])
print(f"pick on tune: {pick}; held-out right {r}/{n} ({100*r/n:.1f}%) vs 0.30: {score('0.30', slugs[half:])[0]}/{n}")
