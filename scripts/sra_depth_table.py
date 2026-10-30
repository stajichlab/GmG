#!/usr/bin/env python3
"""Depth curves for the real-reads validation (#42): accuracy against coverage, paired against single, per exclusion mode.

usage: sra_depth_table.py SRA_RUNS.tsv TRUTH.tsv OUT.md MODE:KIND:N=RESULTS.tsv [MODE:KIND:N=RESULTS.tsv ...]
  MODE  exclusion of the run's own genome (`self`), its species (`species`) or its genus (`genus`) from the reference; KIND `pair` (R1+R2) or `single` (R1);
  N     read pairs used (the first N of each file); RESULTS.tsv: `gmg run` of that setting.
Coverage of an input = bases read / genome size = N x read length x (2 for pairs, 1 for single) / the assembly size of the run's own genome.
Per mode and kind, inputs are pooled over depths and binned by that coverage (<0.1, 0.1-0.25, 0.25-0.5, 0.5-1, 1-2, 2-5, 5-10, >=10). For each bin: n, share with a
genus call, precision of those calls, recall (correct genus calls / n), and the same for species (where the truth species is not hidden). The minimum depth for a
rank is the lowest bin in which recall reaches 0.8 (and stays above it in every deeper bin). Prints Markdown to OUT.md.
"""
import csv, sys, collections

runs = {r["run"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
truth = {r["id"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
BINS = [(0, 0.1, "<0.1x"), (0.1, 0.25, "0.1-0.25x"), (0.25, 0.5, "0.25-0.5x"), (0.5, 1, "0.5-1x"), (1, 2, "1-2x"), (2, 5, "2-5x"), (5, 10, "5-10x"), (10, 1e9, ">=10x")]
bin_of = lambda c: next(lab for lo, hi, lab in BINS if lo <= c < hi)
data = collections.defaultdict(list)   # (mode, kind) -> [(coverage, result row, truth row)]
for spec in sys.argv[4:]:
    key, path = spec.split("=", 1)
    mode, kind, n = key.split(":"); n = int(n)
    for r in csv.DictReader(open(path), delimiter="\t"):
        run = r["query_id"]
        if run not in runs or run not in truth: continue
        gs = int(runs[run]["genome_size"] or 0)
        if gs <= 0: continue
        cov = n * int(runs[run]["read_len"]) * (2 if kind == "pair" else 1) / gs
        data[(mode, kind)].append((cov, n, r, truth[run]))
out = []
for (mode, kind), rows in sorted(data.items()):
    out.append(f"\n### {mode} hidden, {'R1 + R2' if kind == 'pair' else 'R1 only'}  ({len({t['id'] for _, _, _, t in rows})} runs, {len(rows)} inputs)\n")
    out.append("| coverage | n | genus called | genus precision | genus recall | species called | species precision | species recall | phylum precision (called) |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    rec = {}
    for lo, hi, lab in BINS:
        g = [(r, t) for c, _, r, t in rows if lo <= c < hi]
        if not g: continue
        gc = [(r, t) for r, t in g if r["call_genus"]]
        gok = sum(r["call_genus"] == t["genus"] for r, t in gc)
        sc = [(r, t) for r, t in g if r["call_species"]]
        sok = sum(r["call_species"] == t["species"] for r, t in sc)
        pc = [(r, t) for r, t in g if r["call_phylum"]]
        pok = sum(r["call_phylum"] == t["phylum"] for r, t in pc)
        rec[lab] = (gok / len(g), sok / len(g))
        out.append(f"| {lab} | {len(g)} | {len(gc) / len(g):.2f} | {(gok / len(gc)) if gc else float('nan'):.3f} ({len(gc)}) | {gok / len(g):.2f} | {len(sc) / len(g):.2f} | "
                   f"{(sok / len(sc)) if sc else float('nan'):.3f} ({len(sc)}) | {sok / len(g):.2f} | {(pok / len(pc)) if pc else float('nan'):.3f} ({len(pc)}) |")
    labs = [lab for _, _, lab in BINS if lab in rec]
    for i, (rank, idx) in enumerate((("genus", 0), ("species", 1))):
        mins = next((lab for j, lab in enumerate(labs) if all(rec[l][idx] >= 0.8 for l in labs[j:])), None)
        out.append(f"\nMinimum depth for {rank} (recall >= 0.8 from there on): " + (mins or "not reached"))
open(sys.argv[3], "w").write("\n".join(out) + "\n")
print("\n".join(out))
