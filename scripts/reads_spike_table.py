#!/usr/bin/env python3
"""Retention table for the contaminated-reads test (#45): per donor and share of foreign read pairs, how many inputs keep the call they have without contamination.

usage: reads_spike_table.py TRUTH.tsv LABEL=RESULTS.tsv [LABEL=RESULTS.tsv ...]
TRUTH.tsv: id, phylum..species, host, donor, frac (reads_spike_plan.py). For each run (e.g. default, abund1) and each (donor, frac):
n inputs; genus called and correct among all (recall) and among the hosts whose genus call is correct without contamination ("kept");
species correct (hosts with a species truth, i.e. Coccidioides) among the hosts whose species call is correct alone; median top-hit ANI;
flags not_fungal / distant_non_fungal_hit; wrong genus calls (a call that differs from the truth).
"""
import csv, statistics, sys, collections

truth = {r["id"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
for spec in sys.argv[2:]:
    label, path = spec.split("=", 1)
    res = {r["query_id"]: r for r in csv.DictReader(open(path), delimiter="\t")}
    base = {t["host"]: res.get(i) for i, t in truth.items() if t["donor"] == "none" and res.get(i)}
    ok_g = {h for h, r in base.items() if r["call_genus"] and r["call_genus"] == truth[f"{h}__none_0"]["genus"]}
    ok_s = {h for h, r in base.items() if r["call_species"] and truth[f"{h}__none_0"]["species"] and r["call_species"] == truth[f"{h}__none_0"]["species"]}
    print(f"\n### {label}: hosts {len(base)}, genus correct alone {len(ok_g)}, species correct alone {len(ok_s)}\n")
    print("| donor | foreign pairs | n | genus kept | species kept | wrong genus | median ANI | not_fungal | distant hit |")
    print("|---|---|---|---|---|---|---|---|---|")
    groups = collections.defaultdict(list)
    for i, t in truth.items():
        if i in res: groups[(t["donor"], int(t["frac"]))].append((res[i], t))
    for (d, f), L in sorted(groups.items()):
        gk = sum(1 for r, t in L if t["host"] in ok_g and r["call_genus"] == t["genus"])
        sk = sum(1 for r, t in L if t["host"] in ok_s and r["call_species"] == t["species"])
        wrong = sum(1 for r, t in L if r["call_genus"] and r["call_genus"] != t["genus"] and t["host"] in ok_g)
        anis = [float(r["ani"]) for r, t in L if r["ani"]]
        nf = sum("not_fungal" in r["qc_flags"] for r, t in L)
        dh = sum("distant_non_fungal_hit" in r["qc_flags"] for r, t in L)
        print(f"| {d} | {f}% | {len(L)} | {gk}/{len(ok_g)} | {sk}/{len(ok_s)} | {wrong} | {statistics.median(anis) if anis else float('nan'):.4f} | {nf} | {dh} |")
