#!/usr/bin/env python3
"""P2: per-run table and summary statistics from the phase-1 pilot's raw files in results/.

usage: summarize.py   (run in bench/P2_phase1_pilot; reads results/*.tsv, results/downloads_phase1.csv, manifest.csv)
Prints a markdown table and the medians quoted in README.md."""
import csv
import statistics as st

man = {r["run_accession"]: r for r in csv.DictReader(open("manifest.csv"))}
ok = {r["run"]: r for r in csv.DictReader(open("results/results.tsv"), delimiter="\t") if r["outcome"] == "ok"}
qc = {r["id"]: r for r in csv.DictReader(open("results/qc.tsv"), delimiter="\t")}
hook = {r[0]: (int(r[1]), int(r[2])) for r in csv.reader(open("results/hook_times.tsv"), delimiter="\t")}
rss, wall = {}, {}
for line in open("results/sacct.tsv"):
    f = line.strip().split("|")
    if len(f) >= 4 and f[0].endswith(".batch") and f[1] == "COMPLETED":
        key = f[0].split(".")[0].split("_")[-1]
        wall[key], rss[key] = f[2], int(f[3].rstrip("K")) / 1e6
# task index of each run in the manifest; the rerun (job 31992298, task 2) replaces the failed first try
idx = {r["run_accession"]: str(i) for i, r in enumerate(csv.DictReader(open("manifest.csv")))}
dl = {}
for r in csv.reader(open("results/downloads_phase1.csv")):
    if len(r) >= 13 and r[10] == "ok" and r[1] in man:
        dl[r[1]] = (r[4], float(r[8]) if r[8] else None, int(r[6]) if r[6] else None)


def mins(t):
    h, m, s = t.split(":")
    return int(h) * 60 + int(m) + int(s) / 60


print("| Run | Depth | Route | Kept bases | Reached | Transfer | Rate MB/s | Task wall | Peak RSS GB | Sketch s | Sketch MB | Est. size Mb |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
aws_rates, walls, rsss, sk_s, sk_mb = [], [], [], [], []
for run in man:
    r, d, h, q = ok[run], dl[run], hook[run], qc[run]
    w = wall.get(idx[run])
    rate = d[1]
    print(f"| {run} | {float(man[run]['est_depth_x']):.0f}x | {d[0]} | {int(r['bases_kept']) / 1e6:.0f} Mb | "
          f"{float(r['coverage_reached']):.3f}x | {d[2] / 1e6:.0f} MB | {rate:.1f} | {w} | {rss.get(idx[run], float('nan')):.2f} | "
          f"{h[0]} | {h[1] / 1e6:.1f} | {float(q['est_size_mb']):.1f} |")
    if d[0] == "aws_odp":
        aws_rates.append(rate)
    walls.append(mins(w))
    rsss.append(rss[idx[run]])
    sk_s.append(h[0])
    sk_mb.append(h[1] / 1e6)
print()
print(f"AWS rate per task: n={len(aws_rates)}, min {min(aws_rates):.1f}, median {st.median(aws_rates):.1f}, max {max(aws_rates):.1f} MB/s")
print(f"Task wall (all 10): min {min(walls):.1f}, median {st.median(walls):.1f}, max {max(walls):.1f} min")
print(f"Peak RSS (all 10): min {min(rsss):.2f}, median {st.median(rsss):.2f}, max {max(rsss):.2f} GB")
print(f"Sketch step: {min(sk_s)}-{max(sk_s)} s; sketch files per run: min {min(sk_mb):.1f}, median {st.median(sk_mb):.1f}, max {max(sk_mb):.1f} MB")
print(f"Coverage reached: {min(float(ok[r]['coverage_reached']) for r in ok):.3f} to {max(float(ok[r]['coverage_reached']) for r in ok):.3f}x")
