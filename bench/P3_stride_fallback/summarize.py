#!/usr/bin/env python3
"""P3: table of the stride and windows-with-fallback test arrays, from the raw files beside this script.
usage: summarize.py    (reads stride/ and windows/, sacct.tsv and ../P2_phase1_pilot/results/results.tsv)"""
import csv

p2 = {r["run"]: r for r in csv.DictReader(open("../P2_phase1_pilot/results/results.tsv"), delimiter="\t") if r["outcome"] == "ok"}
rss = {}
for line in open("sacct.tsv"):
    f = line.strip().split("|")
    if len(f) >= 4:
        rss[f[0].split(".")[0]] = (f[2], int(f[3].rstrip("K")) / 1e6)
man = [r["run_accession"] for r in csv.DictReader(open("manifest_test.csv"))]
print("| Run | Mode | Sampling | Fell back | Check on the first sample | Check on the final sample | Result s | Phase 1 windows s | Task wall | Peak RSS GB |")
print("|---|---|---|---|---|---|---|---|---|---|")
for mode, job in (("stride", "31994474"), ("windows", "31994475")):
    res = {r["run"]: r for r in csv.DictReader(open(f"{mode}/results.tsv"), delimiter="\t")}
    checks = {}
    for r in csv.DictReader(open(f"{mode}/sampling_check.tsv"), delimiter="\t"):
        checks.setdefault(r["id"], []).append(r)
    for i, run in enumerate(man):
        r, c = res[run], checks[run]
        fmt = lambda x: f"{float(x['est_size_mb']):.1f} Mb, {float(x['kmer_coverage']):.0f}x, {x['verdict']}"
        first, final = fmt(c[0]), fmt(c[-1])
        wall, mem = rss[f"{job}_{i}"]
        print(f"| {run} | {mode} | {r['sampling']} | {r['fallback_from'] or 'no'} | {first} | {final if len(c) > 1 else 'same'} | "
              f"{float(r['seconds']):.0f} | {float(p2[run]['seconds']):.0f} | {wall} | {mem:.2f} |")
