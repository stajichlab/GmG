#!/usr/bin/env python3
"""Tables for reports/validation/reads_calibration.md (#32) from the outputs of 32_reads_calibration.sbatch and 33_coverage_sim.sbatch.

usage: reads_calibration_summary.py CAL_DIR ASSEMBLY_RESULTS.tsv > tables.md
  CAL_DIR              $GR/reads/calibration (holds <ds>/n<N>_a<cutoff>/{results.tsv,results.jsonl,v/validation.tsv}, sim/<genome>/..., *.truth.split.tsv)
  ASSEMBLY_RESULTS.tsv `gmg run` of the ZyGoLife assemblies (query_id = strain, total_bases = assembly size): the genome size of a read set whose
                       id ends in _<strain>
Prints Markdown: (A) fixed 1 / fixed 2 / auto per collection, depth and half (cal = calibration, eval = evaluation): genus and species called /
precision; the cutoffs auto chose; (B) accuracy of coverage_est on simulated reads of known coverage and on real reads of known depth.
"""
import csv, glob, json, os, re, statistics, sys

cal_dir, asm_path = sys.argv[1], sys.argv[2]
RANK = {"zygo": ["phylum", "genus"], "cocci": ["genus", "species"], "rna": ["genus", "species"]}
DEPTH = {"zygo": [15000, 75000, 300000, 1500000], "cocci": [100000, 500000, 2000000], "rna": [100000, 500000, 2000000]}
NAME = {"zygo": "ZyGoLife DNA reads (genus truth = lab label)", "cocci": "Coccidioides DNA reads (species truth)", "rna": "BFD RNA-seq (species truth)"}
CUTS = ["1", "2", "auto"]


def vt(path):
    d = {}
    if not os.path.exists(path):
        return d
    for r in csv.DictReader(open(path), delimiter="\t"):
        d[(r["group"], r["rank"])] = r
    return d


def cell(d, group, rank):
    r = d.get((group, rank))
    if not r:
        return "-"
    called = float(r["called"])
    n = int(r["n"])
    p = r["precision"]
    return f"{called:.2f} / {float(p):.3f} ({round(called * n)}/{n})" if p not in ("", "nan") else f"{called:.2f} / - (0/{n})"


out = ["## A. Abundance cutoff: fixed 1, fixed 2, auto\n",
       "Cells: share of inputs with a call / precision of the calls (calls/inputs). `cal` = calibration half, `eval` = evaluation half (SHA-1 of the sample id).\n"]
for ds in ("zygo", "cocci", "rna"):
    ranks = RANK[ds]
    out.append(f"\n**{NAME[ds]}**\n")
    head = "| reads per file | cutoff | " + " | ".join(f"cal {r}" for r in ranks) + " | " + " | ".join(f"eval {r}" for r in ranks) + " |"
    out += [head, "|" + "---|" * (2 + 2 * len(ranks))]
    for n in DEPTH[ds]:
        for c in CUTS:
            base = f"{cal_dir}/{ds}/n{n}_a{c}"
            d = vt(f"{base}/v/validation.tsv")
            if not d:
                continue
            cells = [cell(d, "split=cal", r) for r in ranks] + [cell(d, "split=eval", r) for r in ranks]
            out.append(f"| {n} | {c} | " + " | ".join(cells) + " |")
    # what auto chose
    chosen = {}
    for n in DEPTH[ds]:
        p = f"{cal_dir}/{ds}/n{n}_aauto/results.jsonl"
        if not os.path.exists(p):
            continue
        cnt = {}
        for line in open(p):
            try:
                v = json.loads(line)["qc"]["abundance_hist"]["min_abund"]
            except Exception:
                continue
            cnt[v] = cnt.get(v, 0) + 1
        chosen[n] = ", ".join(f"{k}: {v}" for k, v in sorted(cnt.items()))
    if chosen:
        out.append("\nCutoffs chosen by `auto` (cutoff: number of read sets): " + "; ".join(f"{n} reads: {s}" for n, s in chosen.items()) + "\n")

# ---------------------------------------------------------------- B. coverage accuracy
out.append("\n## B. Accuracy of `coverage_est`\n")
out.append("Ratio = estimated / true coverage; `within 30%` = share of read sets with ratio in 0.7-1.3.\n")


def stats(label, ratios, extra=""):
    if not ratios:
        return f"| {label} | 0 | - | - | - | - |{extra}"
    ratios = sorted(ratios)
    q = lambda p: ratios[min(len(ratios) - 1, int(p * len(ratios)))]
    within = sum(0.7 <= x <= 1.3 for x in ratios) / len(ratios)
    return f"| {label} | {len(ratios)} | {statistics.median(ratios):.2f} | {q(0.1):.2f}-{q(0.9):.2f} | {within:.2f} |{extra}"


hdr = ["| set | n | median ratio | 10-90% range | within 30% |", "|---|---|---|---|---|"]
# simulated reads of known coverage
out.append("\n**Simulated reads (production scaled 10000), default cutoff 2**\n")
rows = hdr[:]
sim = {}
for tdir in sorted(glob.glob(f"{cal_dir}/sim/*/")):
    tr = os.path.join(tdir, "truth.tsv")
    rs = os.path.join(tdir, "results_a2.tsv")
    if not (os.path.exists(tr) and os.path.exists(rs)):
        continue
    truth = {r["id"]: r for r in csv.DictReader(open(tr), delimiter="\t")}
    for r in csv.DictReader(open(rs), delimiter="\t"):
        t = truth.get(r["query_id"])
        if not t:
            continue
        est = r.get("coverage_est", "")
        key = (float(t["coverage"]), t["error"])
        sim.setdefault(key, []).append(float(est) / float(t["coverage"]) if est else None)
none_cov = {}
for (cov, err), vals in sorted(sim.items()):
    ok = [v for v in vals if v is not None]
    none_cov[(cov, err)] = len(vals) - len(ok)
for cov in sorted({k[0] for k in sim}):
    ratios = [v for (c, e), vs in sim.items() if c == cov for v in vs if v is not None]
    n_none = sum(1 for (c, e), vs in sim.items() if c == cov for v in vs if v is None)
    rows.append(stats(f"{cov:g}x (3 genomes x 3 error rates)", ratios) + (f" no estimate: {n_none}" if n_none else ""))
for err in sorted({k[1] for k in sim}):
    ratios = [v for (c, e), vs in sim.items() if e == err for v in vs if v is not None]
    rows.append(stats(f"error {float(err) * 100:g}% (all coverages)", ratios))
out += rows

# real reads: ZyGoLife strains with an assembly, Coccidioides at a nominal genome size
asm = {}
if os.path.exists(asm_path):
    for r in csv.DictReader(open(asm_path), delimiter="\t"):
        try:
            asm[r["query_id"].replace(".sorted", "")] = float(r["total_bases"])
        except ValueError:
            pass
out.append("\n**Real reads, default cutoff 2: coverage = sketched bases / genome size**\n")
rows = hdr[:]
for ds, size_of in (("zygo", lambda rid: next((asm[a] for a in asm if rid.endswith("_" + a)), None)), ("cocci", lambda rid: 28.5e6)):
    for n in DEPTH[ds]:
        p = f"{cal_dir}/{ds}/n{n}_a2/results.tsv"
        if not os.path.exists(p):
            continue
        ratios = []
        for r in csv.DictReader(open(p), delimiter="\t"):
            size = size_of(r["query_id"])
            if size and r.get("coverage_est") and r.get("total_bases"):
                true = float(r["total_bases"]) / size
                if true > 0:
                    ratios.append(float(r["coverage_est"]) / true)
        rows.append(stats(f"{ds}, first {n} reads per file (true coverage from {'its assembly' if ds == 'zygo' else 'a nominal 28.5 Mb genome'})", ratios))
out += rows
print("\n".join(out))
