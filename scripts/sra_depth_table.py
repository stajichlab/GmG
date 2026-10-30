#!/usr/bin/env python3
"""Depth curves for the real-reads validation (#42): accuracy against depth, paired against single, per exclusion mode, per clade and lab.

usage: sra_depth_table.py SRA_RUNS.tsv TRUTH.tsv OUT.md MODE:KIND:N=RESULTS.tsv [MODE:KIND:N=RESULTS.tsv ...]
  MODE  exclusion of the run's own genome (`self`), its species (`species`) or its genus (`genus`) from the reference; KIND `pair` (R1+R2) or `single` (R1);
  N     read pairs used (the first N of each file); RESULTS.tsv: `gmg run` of that setting.
SRA_RUNS.tsv has, per run, the number of OTHER reference genomes of the same species and genus (species_others_in_ref, genus_others_in_ref). A rank is *attainable* for a
run when the truth can still be found in the reference after the exclusion: with the run's own genome hidden the genus is attainable if another genome of the genus
exists and the species if another genome of the species exists; with the species hidden the genus is attainable if another species of the genus exists and the species never
is; with the genus hidden neither is. **Recall is counted over attainable inputs only** (a run whose genus is nowhere else in the reference cannot be right, whatever the depth);
precision is over the calls made. Tables: by read pairs N (the same for every run), by coverage (N x read length x 2 or 1 / genome size, pooled over runs, so it mixes
genome sizes), by phylum and by submitting center at the deepest N. The minimum depth for a rank is the lowest N (and its median coverage) from which recall over the attainable
inputs stays >= 0.8, or "not reached".
"""
import csv, statistics, sys, collections

runs = {r["run"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
truth = {r["id"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
BINS = [(0, 0.1, "<0.1x"), (0.1, 0.25, "0.1-0.25x"), (0.25, 0.5, "0.25-0.5x"), (0.5, 1, "0.5-1x"), (1, 2, "1-2x"), (2, 5, "2-5x"), (5, 10, "5-10x"), (10, 1e9, ">=10x")]


def attainable(mode, run):
    so, go = int(run["species_others_in_ref"]), int(run["genus_others_in_ref"])
    if mode == "self": return go > 0, so > 0
    if mode == "species": return go - so > 0, False
    return False, False


data = collections.defaultdict(list)   # (mode, kind) -> [(coverage, N, result, truth, run)]
for spec in sys.argv[4:]:
    key, path = spec.split("=", 1)
    mode, kind, n = key.split(":"); n = int(n)
    for r in csv.DictReader(open(path), delimiter="\t"):
        rid = r["query_id"]
        if rid not in runs or rid not in truth or int(runs[rid]["genome_size"] or 0) <= 0: continue
        cov = n * int(runs[rid]["read_len"]) * (2 if kind == "pair" else 1) / int(runs[rid]["genome_size"])
        data[(mode, kind)].append((cov, n, r, truth[rid], runs[rid]))


def stats(rows, mode):
    """(n, n attainable genus, genus called, genus right, genus recall over attainable, same for species, phylum calls, phylum right)"""
    ga = [x for x in rows if attainable(mode, x[4])[0]]
    sa = [x for x in rows if attainable(mode, x[4])[1]]
    gc = [x for x in rows if x[2]["call_genus"]]
    gr = sum(x[2]["call_genus"] == x[3]["genus"] for x in gc)
    gra = sum(x[2]["call_genus"] == x[3]["genus"] for x in ga)
    sc = [x for x in rows if x[2]["call_species"]]
    sr = sum(x[2]["call_species"] == x[3]["species"] for x in sc)
    sra = sum(x[2]["call_species"] == x[3]["species"] for x in sa)
    gta = sum(x[2]["call_genus"] == x[3]["genus"] or x[2].get("tentative_genus", "") == x[3]["genus"] for x in ga)
    pc = [x for x in rows if x[2]["call_phylum"]]
    pr = sum(x[2]["call_phylum"] == x[3]["phylum"] for x in pc)
    return dict(n=len(rows), ga=len(ga), gc=len(gc), gr=gr, grec=gra / len(ga) if ga else float("nan"), sa=len(sa), sc=len(sc), sr=sr,
                srec=sra / len(sa) if sa else float("nan"), pc=len(pc), pr=pr,
                gt=gta / len(ga) if ga else float("nan"))


def row(label, s):
    f = lambda a, b: f"{a / b:.3f} ({b})" if b else "- (0)"
    return (f"| {label} | {s['n']} | {s['ga']} | {s['gc'] / s['n']:.2f} | {f(s['gr'], s['gc'])} | {s['grec']:.2f} | {s['sa']} | {s['sc'] / s['n']:.2f} | "
            f"{f(s['sr'], s['sc'])} | {s['srec']:.2f} | {f(s['pr'], s['pc'])} | {s['gt']:.2f} |")


HEAD = ("| {k} | n | genus attainable | genus called | genus precision | genus recall (attainable) | species attainable | species called | species precision | "
        "species recall (attainable) | phylum precision | genus or tentative genus right (attainable) |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
out = []
mins = {}
for (mode, kind), rows in sorted(data.items()):
    runs_n = len({x[3]["id"] for x in rows})
    out.append(f"\n### {mode} hidden, {'R1 + R2' if kind == 'pair' else 'R1 only'} ({runs_n} runs)\n\n**By read pairs N**\n\n" + HEAD.format(k="N"))
    byn = collections.defaultdict(list)
    for x in rows: byn[x[1]].append(x)
    rec = {}
    for n in sorted(byn):
        s = stats(byn[n], mode); rec[n] = s
        med = statistics.median(x[0] for x in byn[n])
        out.append(row(f"{n:,} (median {med:.2f}x)", s))
    for rank, key in (("genus", "grec"), ("species", "srec")):
        ns = sorted(rec)
        first = next((n for i, n in enumerate(ns) if rec[n][key] == rec[n][key] and all(rec[m][key] >= 0.8 for m in ns[i:])), None)
        mins[(mode, kind, rank)] = (first, statistics.median(x[0] for x in byn[first]) if first else None)
        best = max((rec[n][key] for n in ns if rec[n][key] == rec[n][key]), default=float("nan"))
        out.append(f"\nMinimum depth for {rank} (recall over attainable inputs >= 0.8 from there on): "
                   + (f"{first:,} pairs (median {mins[(mode, kind, rank)][1]:.2f}x)" if first else f"not reached (best {best:.2f})"))
    out.append("\n**By coverage (pooled over runs and N; coverage and genome size/clade are confounded)**\n\n" + HEAD.format(k="coverage"))
    for lo, hi, lab in BINS:
        g = [x for x in rows if lo <= x[0] < hi]
        if g: out.append(row(lab, stats(g, mode)))
    if kind == "pair" and mode in ("self", "species"):
        top = max(byn)
        deep = byn[top]
        out.append(f"\n**By phylum, {top:,} pairs**\n\n" + HEAD.format(k="phylum"))
        for ph in sorted({x[3]["phylum"] for x in deep}, key=lambda p: -sum(x[3]["phylum"] == p for x in deep)):
            out.append(row(ph, stats([x for x in deep if x[3]["phylum"] == ph], mode)))
        out.append(f"\n**By submitting center (centers with >= 3 runs), {top:,} pairs**\n\n" + HEAD.format(k="center"))
        cen = collections.Counter(x[3]["center"] for x in deep)
        for c, k in cen.most_common():
            if k >= 3: out.append(row(c, stats([x for x in deep if x[3]["center"] == c], mode)))
out.append("\n### Minimum depth, summary\n\n| mode | reads | genus | species |\n|---|---|---|---|")
for (mode, kind) in sorted({(m, k) for m, k, _ in mins}):
    cell = lambda r: (f"{mins[(mode, kind, r)][0]:,} pairs ({mins[(mode, kind, r)][1]:.2f}x)" if mins[(mode, kind, r)][0] else "not reached")
    out.append(f"| {mode} hidden | {'R1 + R2' if kind == 'pair' else 'R1 only'} | {cell('genus')} | {cell('species')} |")
open(sys.argv[3], "w").write("\n".join(out) + "\n")
print("\n".join(out))
