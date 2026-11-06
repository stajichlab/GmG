#!/usr/bin/env python3
"""Run a fixed popsketch on the finished sketch sets and compare with the earlier run of the same sketches (<slug>_out).
usage (from popsra): NEW=path/to/new/popsketch 09_validate_estimator.py SLUG [SLUG ...]   (new results go to <slug>_fix)"""
import csv, json, os, statistics, subprocess, sys

NEW = os.environ["NEW"]
def load(d):
    q = list(csv.DictReader(open(f"{d}/qc.tsv"), delimiter="\t"))
    return q, json.load(open(f"{d}/popsketch.json"))
for s in sys.argv[1:]:
    subprocess.run([NEW, "run", "--from-sketches", f"{s}/sketches", "-o", f"{s}_fix", "--threads", "8", "-q"], check=True)
    (qo, jo), (qn, jn) = load(f"{s}_out"), load(f"{s}_fix")
    def summ(q, j):
        basis = {}
        for r in q: basis[r["coverage_basis"] or "none"] = basis.get(r["coverage_basis"] or "none", 0) + 1
        size = [float(r["est_size_mb"]) for r in q if r["est_size_mb"]]
        g = j["clusters"]["gap"]
        fl = {k: len(v) for k, v in j["flags"].items()}
        return f'in_analysis {j["samples"]["in_analysis"]}/{j["samples"]["total"]}  basis {basis}  size Mb median {statistics.median(size):.1f} (p10 {sorted(size)[len(size)//10]:.1f}, p90 {sorted(size)[9*len(size)//10]:.1f})  {g["verdict"]} jump {g["gap_over_floor"]:.2f}  floor {j["noise_floor"]["value"]:.4f}  flags {fl}'
    print(f"{s}\n  old: {summ(qo, jo)}\n  new: {summ(qn, jn)}")
    ch = [(a["id"], a["coverage_basis"], b["coverage_basis"], a["est_size_mb"], b["est_size_mb"]) for a, b in zip(qo, qn) if a["id"] == b["id"] and (a["coverage_basis"] != b["coverage_basis"])]
    print(f"  {len(ch)} samples changed basis" + (": " + "; ".join(f"{i} {x}->{y} ({sa}->{sb} Mb)" for i, x, y, sa, sb in ch[:6]) if ch else ""))
