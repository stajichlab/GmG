#!/usr/bin/env python3
"""Flag rate and error estimate of `high_error_reads` on the re-simulated reads of scripts/hpcc/39_long_reads_flag_check.sbatch (#44).

usage: long_reads_flag_summary.py FLAG_DIR OUT.tsv   (FLAG_DIR/<genome>/out/results.jsonl; ids are e<error>_d<depth>)
Writes one row per (error, depth): n, flagged share, median error_est, share with no estimate; prints a markdown table of the flag rate.
"""
import glob, json, sys, statistics, collections

rows = collections.defaultdict(list)
for p in glob.glob(f"{sys.argv[1]}/*/out/results.jsonl"):
    for line in open(p):
        r = json.loads(line)
        e, d = r["id"][1:].split("_d")
        q = r["qc"] or {}
        rows[(float(e) * 100, float(d))].append(("high_error_reads" in r["qc_flags"], q.get("error_est")))
with open(sys.argv[2], "w") as f:
    f.write("error_pct\tdepth\tn\tflag_share\tmedian_error_est\tno_estimate_share\n")
    for (e, d), v in sorted(rows.items()):
        est = [x[1] for x in v if x[1] is not None]
        f.write(f"{e:g}\t{d:g}\t{len(v)}\t{sum(x[0] for x in v)/len(v):.3f}\t{statistics.median(est) if est else ''}\t{1-len(est)/len(v):.3f}\n")
deps = sorted({d for _, d in rows})
print("| error \\ depth | " + " | ".join(f"{d:g}x" for d in deps) + " |\n|---|" + "---|" * len(deps))
for e in sorted({e for e, _ in rows}):
    cells = []
    for d in deps:
        v = rows[(e, d)]
        est = [x[1] for x in v if x[1] is not None]
        cells.append(f"{sum(x[0] for x in v)}/{len(v)}" + (f" ({statistics.median(est)*100:.1f}%)" if est else ""))
    print(f"| {e:g}% | " + " | ".join(cells) + " |")
