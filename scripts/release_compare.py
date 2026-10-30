#!/usr/bin/env python3
"""What changes between two index releases, and does the calibration still hold (#48).

usage: release_compare.py OUT.md --order V1,V2,... --set NAME=truth.tsv|- [--excluded FILE] [--decoys meta_decoys.tsv]
                          --results VER:NAME=results.tsv ...
  --results   the results.tsv of one test set classified against one release (same binary, same thresholds)
  --set       the test set's truth table (id + phylum..species columns; `-` = no truth: only the change counts are made)
  --excluded  ids of reference genomes removed from the reference in some release (reports/label_audit/reference_exclusions.txt)
  --decoys    meta_decoys.tsv (ids of the non-fungal decoy genomes)
For every consecutive pair of releases and every set: queries whose status or call at any rank changed, with the cause of each change in
this order: the old top hit was excluded from the new release; a decoy is now (or no longer) the top hit; the same top hit with another
label (taxonomy update); another nearest genome; other (thresholds, flags). With a truth table: precision and called share per rank per
release, and the RE-CALIBRATION TRIGGER (docs/recalibration.md): exit code 3 and a line `TRIGGER` when a rule fires.
"""
import argparse, collections, csv, math, sys

RANKS = ["phylum", "class", "order", "family", "genus", "species"]
ap = argparse.ArgumentParser()
ap.add_argument("out")
ap.add_argument("--order", required=True)
ap.add_argument("--set", action="append", default=[])
ap.add_argument("--results", action="append", default=[])
ap.add_argument("--excluded")
ap.add_argument("--decoys")
a = ap.parse_args()
order = a.order.split(",")
sets = {}
for s in a.set:
    n, t = s.split("=", 1)
    sets[n] = None if t == "-" else {r["id"]: r for r in csv.DictReader(open(t), delimiter="\t")}
res = collections.defaultdict(dict)  # res[set][version] = {query_id: row}
for r in a.results:
    k, p = r.split("=", 1)
    v, n = k.split(":", 1)
    res[n][v] = {x["query_id"]: x for x in csv.DictReader(open(p), delimiter="\t")}
excl = {l.strip() for l in open(a.excluded)} if a.excluded else set()
decoys = {r["id"] for r in csv.DictReader(open(a.decoys), delimiter="\t")} if a.decoys else set()


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def prec(rows, truth, rank):
    t = [(r, truth[q]) for q, r in rows.items() if q in truth and r["status"] != "error" and truth[q].get(rank)]
    called = [(r, g) for r, g in t if r[f"call_{rank}"]]
    right = sum(1 for r, g in called if r[f"call_{rank}"] == g[rank])
    return len(t), len(called), right


def top_id(r):
    return r.get("top_hit", "")


def cause(o, n):
    th_o, th_n = top_id(o), top_id(n)
    if th_o and th_o in excl:
        return "old top hit excluded from the reference"
    if (th_n in decoys) != (th_o in decoys):
        return "a decoy genome is (no longer) the top hit"
    if th_o and th_o == th_n:
        return "same top hit, new label (taxonomy update)"
    if th_o != th_n:
        return "another nearest genome"
    return "other (thresholds, flags)"


md = ["# Release stability (#48)\n", f"Releases compared: {', '.join(order)}. Same binary and thresholds for all runs (`known`).\n"]
trigger = []
notes = []
for name in sorted(res):
    vs = [v for v in order if v in res[name]]
    md.append(f"## {name}\n")
    # changes between consecutive releases
    md.append("| from | to | inputs | status changed | any call changed | " + " | ".join(RANKS) + " | `ok` in both, genus or species call changed |")
    md.append("|---|---|---|---|---|" + "---|" * len(RANKS) + "---|")
    causes = {}
    for v0, v1 in zip(vs, vs[1:]):
        A, B = res[name][v0], res[name][v1]
        ids = [q for q in A if q in B]
        st = sum(1 for q in ids if A[q]["status"] != B[q]["status"])
        per = {r: sum(1 for q in ids if A[q][f"call_{r}"] != B[q][f"call_{r}"]) for r in RANKS}
        anyc = [q for q in ids if A[q]["status"] != B[q]["status"] or any(A[q][f"call_{r}"] != B[q][f"call_{r}"] for r in RANKS)]
        okc = sum(1 for q in ids if A[q]["status"] == "ok" == B[q]["status"] and (A[q]["call_genus"] != B[q]["call_genus"] or A[q]["call_species"] != B[q]["call_species"]))
        md.append(f"| {v0} | {v1} | {len(ids)} | {st} | {len(anyc)} ({100 * len(anyc) / max(1, len(ids)):.1f}%) | " + " | ".join(str(per[r]) for r in RANKS) + f" | {okc} |")
        causes[(v0, v1)] = collections.Counter(cause(A[q], B[q]) for q in anyc)
    md.append("")
    md.append("Causes of the changed inputs (first matching rule):\n")
    allc = sorted({c for cs in causes.values() for c in cs})
    if allc:
        md.append("| cause | " + " | ".join(f"{a_}→{b_}" for a_, b_ in causes) + " |")
        md.append("|---|" + "---|" * len(causes))
        for c in allc:
            md.append(f"| {c} | " + " | ".join(str(causes[k].get(c, 0)) for k in causes) + " |")
        md.append("")
    # precision per release
    truth = sets.get(name)
    if truth:
        md.append("Precision of the calls made (called share), against the truth table:\n")
        md.append("| release | " + " | ".join(RANKS) + " |")
        md.append("|---|" + "---|" * len(RANKS))
        stats = {}
        for v in vs:
            cells = []
            for r in RANKS:
                n, c, k = prec(res[name][v], truth, r)
                stats[(v, r)] = (n, c, k)
                lo, hi = wilson(k, c)
                cells.append(f"{k / c:.3f} ({c / n:.2f})" if c else "-")
            md.append(f"| {v} | " + " | ".join(cells) + " |")
        md.append("")
        for v0, v1 in zip(vs, vs[1:]):
            for r in RANKS:
                n0, c0, k0 = stats[(v0, r)]
                n1, c1, k1 = stats[(v1, r)]
                if c0 and c1 and k1 / c1 < wilson(k0, c0)[0]:
                    trigger.append(f"{name} {v0}->{v1} {r}: precision {k1 / c1:.3f} is below the 95% lower bound {wilson(k0, c0)[0]:.3f} of the previous release")
                if r in ("genus", "species") and c0 and abs(c1 / n1 - c0 / n0) > 0.05:
                    trigger.append(f"{name} {v0}->{v1} {r}: called share moved by more than 5 points ({c0 / n0:.3f} -> {c1 / n1:.3f})")
    # confident-call churn
    for v0, v1 in zip(vs, vs[1:]):
        A, B = res[name][v0], res[name][v1]
        okc = [q for q in A if q in B and A[q]["status"] == "ok" and B[q]["status"] == "ok"]
        changed = [q for q in okc if A[q]["call_genus"] != B[q]["call_genus"]]
        explained = ("old top hit excluded from the reference", "same top hit, new label (taxonomy update)")
        ch = sum(1 for q in changed if cause(A[q], B[q]) not in explained)
        notes.append(f"{name} {v0}->{v1}: {len(changed)} of {len(okc)} confident genus calls changed, {len(changed) - ch} of them explained "
                     f"(a deliberately excluded top hit or a label update), {ch} unexplained")
        if okc and ch / len(okc) > 0.02:
            trigger.append(f"{name} {v0}->{v1}: {ch} of {len(okc)} confident genus calls changed for no documented reason ({100 * ch / len(okc):.1f}% > 2%)")
md.append("## Re-calibration trigger\n")
md.append("Churn of confident genus calls (rule 3 counts only the unexplained ones):\n")
md += [f"* {n}" for n in notes]
md.append("")
if trigger:
    md.append("TRIGGER: " + "; ".join(trigger) + "\n")
else:
    md.append("No rule fired (see docs/recalibration.md for the rules).\n")
open(a.out, "w").write("\n".join(md) + "\n")
print("\n".join(md))
for t in trigger:
    print("TRIGGER", t, file=sys.stderr)
sys.exit(3 if trigger else 0)
