#!/usr/bin/env python3
"""P4: the numbers of the phase 1 species run, from the raw files in results/ (run in bench/P4_phase1_species).
usage: summarize.py > results/summary.md"""
import collections
import csv
import io
import json
import statistics as st
import tarfile

T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
man = list(csv.DictReader(open("manifest.csv")))
res = T("results/results.tsv")
side = [json.load(io.TextIOWrapper(t.fileobj if False else tarfile.open("results/sidecars.tgz").extractfile(t), "utf-8"))
        for t in tarfile.open("results/sidecars.tgz").getmembers() if t.name.endswith(".json")]
S = {s["run"]: s for s in side}
chunk_node = {0: "r39", 1: "r39", 2: "r11", 3: "r11", 4: "r11", 5: "r11", 6: "r12", 7: "r12", 8: "r12", 9: "r12"}
fetches_on_node = {"r39": 8, "r11": 16, "r12": 8}  # PAR=4 per job x jobs sharing the node
node_of = {r["run_accession"]: chunk_node[i // 32] for i, r in enumerate(man)}


def q(v, p):
    v = sorted(v)
    return v[min(len(v) - 1, int(p * len(v)))]


def line(name, v, unit="", f="{:.1f}"):
    print(f"| {name} | n = {len(v)} | min {f.format(min(v))} | median {f.format(st.median(v))} | p90 {f.format(q(v, 0.9))} | max {f.format(max(v))} {unit}|")


print("## Outcome")
oc = collections.Counter(r["outcome"] for r in res)
ex = collections.Counter(l.split("\t")[1] for l in open("results/exits.tsv"))
print(f"- {len(res)} result lines for {len(man)} manifest runs: {dict(oc)}; exit codes {dict(ex)}")
print(f"- routes {dict(collections.Counter(s['route'] for s in side))}; sampling {dict(collections.Counter(s['sampling'] for s in side))}; "
      f"verify scope {dict(collections.Counter(s['verify_scope'] for s in side))}; fallback {dict(collections.Counter(bool(s['fallback_from']) for s in side))}")
print(f"- declared concurrency in the sidecars: {dict(collections.Counter(s['concurrency'] for s in side))}")
print(f"- coverage reached: min {min(float(r['coverage_reached']) for r in res):.3f}x, median {st.median(float(r['coverage_reached']) for r in res):.3f}x, "
      f"{sum(1 for r in res if float(r['coverage_reached']) < 24.9)} runs below 24.9x (runs shorter than the target are kept whole)")
tot_b = sum(s["bytes_total"] for s in side) / 1e9
print(f"- transferred {tot_b:.0f} GB for {sum(s['bases_kept'] for s in side) / 1e9:.0f} Gb of bases kept from {sum(int(r['base_count']) for r in man) / 1e9:.0f} Gb of runs")

bpb = [S[m["run_accession"]]["bytes_total"] / int(m["base_count"]) for m in man if S[m["run_accession"]]["route"] == "aws_odp"]
print(f"- `.sra` bytes per base over AWS ({len(bpb)} runs): min {min(bpb):.2f}, median {st.median(bpb):.2f}, max {max(bpb):.2f} (the earlier range, n = 3: 0.28 to 0.65)")
print("\n## Per-fetch cost")
print("| quantity | n | min | median | p90 | max |")
print("|---|---|---|---|---|---|")
line("fetch total seconds (results.tsv)", [float(r["seconds"]) for r in res])
line("transfer seconds (sidecar)", [s["seconds_fetch"] for s in side])
line("conversion seconds (sidecar)", [s["seconds_convert"] for s in side])
hk = {h[0]: int(h[1]) for h in (l.rstrip("\n").split("\t") for l in open("results/hook_times.tsv"))}
resid = [float(r["seconds"]) - S[r["run"]]["seconds_fetch"] - S[r["run"]]["seconds_convert"] - hk[r["run"]] for r in res]
share = [x / float(r["seconds"]) for x, r in zip(resid, res)]
print(f"\nTime the sidecar timers do not name (fetch seconds minus transfer, conversion and sketch; per run): median {st.median(resid):.0f} s, "
      f"{100 * st.median(share):.0f}% of the fetch (p10 {100 * q(share, 0.1):.0f}%, p90 {100 * q(share, 0.9):.0f}%). It includes the plan lookups, "
      "the read counts of the full files, the evenly spaced cut and the checks after it; they were not timed separately, so which of them "
      "dominates is UNMEASURED.")
print(f"\nSum of per-fetch seconds: {sum(float(r['seconds']) for r in res) / 3600:.1f} task-hours (my estimate from the pilot: 16 to 33)")
ch = {}
for r in man:
    ch.setdefault(r["run_accession"], r)
sec_per_mb = [float(r["seconds"]) / (int(m["base_count"]) / 1e6) for r in res for m in man if m["run_accession"] == r["run"]]
line("fetch seconds per Mb of run", sec_per_mb, "(estimate used: 0.19 to 0.40)", "{:.2f}")

print("\n## Transfer rate per fetch (MB/s, transfer only, from the sidecars) by how many fetches shared the node")
print("| node | fetches on the node | n | min | median | max |")
print("|---|---|---|---|---|---|")
for node in ("r12", "r39", "r11"):
    v = [S[r]["mb_per_s"] for r in S if node_of[r] == node and S[r]["route"] == "aws_odp" and S[r]["mb_per_s"]]
    print(f"| {node} | {fetches_on_node[node]} | {len(v)} | {min(v):.1f} | {st.median(v):.1f} | {max(v):.1f} |")
v = [S[r]["mb_per_s"] for r in S if S[r]["route"] == "ena_https"]
print(f"\nENA ({len(v)} runs): " + ", ".join(f"{x:.2f}" for x in v) + " MB/s")

print("\n## Jobs (sacct)")
print("| task | node | elapsed | peak RSS GB (of 20 GiB requested) |")
print("|---|---|---|---|")
for l in open("results/sacct.txt"):
    f = l.strip().split("|")
    if len(f) == 5 and f[0].endswith(".batch") and f[1] == "COMPLETED":
        print(f"| {f[0].split('.')[0]} | {f[4]} | {f[2]} | {int(f[3].rstrip('K')) / 1048576:.1f} GiB |")

print("\n## Sketches and the sampling check")
hook = [l.rstrip("\n").split("\t") for l in open("results/hook_times.tsv")]
line("sketch seconds", [int(h[1]) for h in hook], "", "{:.0f}")
line("sketch files per run (MB)", [int(h[2]) / 1e6 for h in hook])
chk = T("results/sampling_check.tsv")
print(f"\n- `popsketch check`: {dict(collections.Counter(c['verdict'] for c in chk))}; size ratio {min(float(c['size_ratio']) for c in chk):.2f} to "
      f"{max(float(c['size_ratio']) for c in chk):.2f} (median {st.median(float(c['size_ratio']) for c in chk):.2f}); coverage ratio "
      f"{min(float(c['coverage_ratio']) for c in chk if c['coverage_ratio']):.2f} to {max(float(c['coverage_ratio']) for c in chk if c['coverage_ratio']):.2f}")

print("\n## The species analysis (popsketch run --collapse mean --collapse-within study)")
qc = T("results/qc.tsv")
fl = collections.Counter()
for r in qc:
    for f in filter(None, (r["flags"] or "").replace(";", ",").split(",")):
        fl[f.strip()] += 1
print(f"- audit: {dict(collections.Counter(r['status'] for r in qc))}; flags {dict(fl)}")
by = collections.defaultdict(lambda: [0, 0])
for r in qc:
    st_ = next(m["study_accession"] for m in man if m["run_accession"] == r["id"])
    by[st_][1] += 1
    by[st_][0] += bool(r["flags"])
big = sorted(by.items(), key=lambda x: -x[1][1])[:3]
print("- flagged / total by the three largest studies: " + ", ".join(f"{s} {a}/{b}" for s, (a, b) in big))
lowc = [r for r in qc if "low_coverage" in (r["flags"] or "")]
ok = [r for r in qc if not r["flags"]]
print(f"- k-mer coverage: low_coverage samples median {st.median(float(r['kmer_coverage']) for r in lowc):.1f} (n = {len(lowc)}), unflagged median {st.median(float(r['kmer_coverage']) for r in ok):.1f} (n = {len(ok)})")
col = T("results/collapse.tsv")
reps = sum(r["is_representative"] == "true" for r in col)
print(f"- collapse: {len(col)} samples to {reps} representatives (largest cluster {max(int(r['cluster_size']) for r in col)}); noise floor in popsketch.json")
print(f"- label tests: ")
for g in T("results/group_tests.tsv"):
    print(f"  - {g['label_column']}: n {g['n']}, groups {g['groups']}, R2 {float(g['r2']):.3f}, p {g['p']}, strata {g['strata']}")
