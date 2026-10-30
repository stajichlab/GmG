#!/usr/bin/env python3
"""Markdown tables for reports/validation/tool_comparison.md. usage: make_table.py COMPARE_DIR SCORE.tsv > tables.md
Reads score.py's output and the GNU time files of the runs (wall time, peak memory), the reference build times and the index sizes."""
import collections, csv, glob, os, re, sys

cdir, score = sys.argv[1:3]
R = ["phylum", "class", "order", "family", "genus", "species"]
rows = list(csv.DictReader(open(score), delimiter="\t"))
by = {(r["tool"], r["variant"], r["set"], r["mode"], r["group"], r["rank"]): r for r in rows}
TV = [("gmg", "native"), ("skani", "raw"), ("skani", "cutoff"), ("mash", "raw"), ("mash", "cutoff"),
      ("mashscreen", "raw"), ("mashscreen", "cutoff"), ("sourmash", "raw"), ("sourmash", "cutoff")]
NAME = {"gmg": "GmG", "skani": "skani", "mash": "mash dist", "mashscreen": "mash screen", "sourmash": "sourmash gather"}
VAR = {"native": "(own calls)", "raw": "nearest genome, no cutoff", "cutoff": "nearest genome, GmG ANI cutoffs"}


def cell(tool, var, s, mode, group, rank):
    r = by.get((tool, var, s, mode, group, rank))
    if not r:
        return "-"
    p = r["precision"]
    return f"{float(p):.3f} ({float(r['called']):.2f})" if p != "" else f"- ({float(r['called']):.2f})"


def label(t, v):
    return f"{NAME[t]} {VAR[v]}"


out = []
out.append("**Table 1. Post-freeze assemblies (s1, n = 300; independent of the reference). Precision of the calls made (share of inputs called).**\n")
out.append("| tool | " + " | ".join(R) + " |\n|---|" + "---|" * 6)
for t, v in TV:
    if any(k[0] == t and k[1] == v and k[2] == "s1" for k in by):
        out.append(f"| {label(t, v)} | " + " | ".join(cell(t, v, "s1", "none", "all", r) for r in R) + " |")
for title, mode in (("leave-one-out", "loo"), ("leave-species-out", "los"), ("leave-genus-out", "log")):
    out.append(f"\n**Table 2{'abc'[['loo', 'los', 'log'].index(mode)]}. RefSeq truth-set genomes, {title} (s2). Precision (share called).**\n")
    out.append("| tool | " + " | ".join(R) + " |\n|---|" + "---|" * 6)
    for t, v in TV:
        if any(k[0] == t and k[1] == v and k[2] == "s2" and k[3] == mode for k in by):
            out.append(f"| {label(t, v)} | " + " | ".join(cell(t, v, "s2", mode, "all", r) for r in R) + " |")
groups = sorted({k[4] for k in by if k[2] == "s3" and k[4] not in ("all", "negative_controls")})
out.append("\n**Table 3. Stress inputs (s3): genus precision (share called) per degradation; the last column is the share of negative controls that got a phylum call (lower is better).**\n")
out.append("| tool | " + " | ".join(groups) + " | negatives called |\n|---|" + "---|" * (len(groups) + 1))
for t, v in TV:
    if any(k[0] == t and k[1] == v and k[2] == "s3" for k in by):
        neg = by.get((t, v, "s3", "none", "negative_controls", "phylum"))
        out.append(f"| {label(t, v)} | " + " | ".join(cell(t, v, "s3", "none", g, "genus") for g in groups) + f" | {float(neg['called']):.2f} |" if neg else " - |")


def tfile(p):
    d = {}
    if os.path.exists(p):
        t = open(p).read()
        m = re.search(r"Elapsed \(wall clock\) time.*?: ([\d:.]+)", t)
        if m:
            parts = [float(x) for x in m.group(1).split(":")]
            d["wall"] = sum(x * 60 ** i for i, x in enumerate(reversed(parts)))
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
        if m:
            d["rss"] = int(m.group(1)) / 1e6
    return d


nq = {s: len(glob.glob(f"{cdir}/q/{s}/*.fa.gz")) for s in ("s1", "s2", "s3")}
out.append("\n**Table 4. Cost (8 CPUs per run; wall seconds per input over the whole run, including loading the reference once; peak memory).**\n")
out.append("| tool | s1 s/input | s2 s/input | s3 s/input | peak RAM GB | reference build (wall) | reference size |\n|---|---|---|---|---|---|---|")
bt = {"skani": "skani_build.time", "mash": "mash_build.time", "mashscreen": "mash_build.time", "sourmash": "sourmash_build.time"}
size = {"skani": f"{cdir}/ref/skani_db", "mash": f"{cdir}/ref/mash_refs.msh", "mashscreen": f"{cdir}/ref/mash_refs.msh", "sourmash": f"{cdir}/ref/sm_refs.zip"}


def du(p):
    if os.path.isdir(p):
        return sum(os.path.getsize(f) for f in glob.glob(p + "/**", recursive=True) if os.path.isfile(f))
    return os.path.getsize(p) if os.path.exists(p) else 0


for t in ["gmg", "skani", "mash", "mashscreen", "sourmash"]:
    per, rss = [], []
    for s in ("s1", "s2", "s3"):
        if t == "gmg" and s == "s2":
            ds = [tfile(f"{cdir}/out/gmg/s2/run.{m}.time") for m in ("loo", "los", "log")]
            w = sum(d.get("wall", 0) for d in ds) / 3
            rss += [d["rss"] for d in ds if "rss" in d]
        else:
            d = tfile(f"{cdir}/out/{t}/{s}/run.time")
            w = d.get("wall", 0)
            if "rss" in d:
                rss.append(d["rss"])
        per.append(f"{w / nq[s]:.2f}" if w else "-")
    b = tfile(f"{cdir}/ref/{bt[t]}") if t in bt else {}
    if t == "gmg":
        bs, sz = "(index build ~10 min on the full BFD; reps index from the release)", f"{os.path.getsize(glob.glob(cdir + '/../release/2026.10.3/bfd-reps-2026.10.3-s20000.gmgi')[0]) / 1e6:.0f} MB"
    else:
        bs = f"{b['wall'] / 60:.0f} min" if b.get("wall") else "-"
        sz = f"{du(size[t]) / 1e6:.0f} MB"
    out.append(f"| {NAME[t]} | {per[0]} | {per[1]} | {per[2]} | {max(rss):.1f} | {bs} | {sz} |" if rss else f"| {NAME[t]} | {' | '.join(per)} | - | {bs} | {sz} |")
print("\n".join(out))
