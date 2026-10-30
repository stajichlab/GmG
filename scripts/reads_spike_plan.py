#!/usr/bin/env python3
"""Hosts and truth for the contaminated-reads test (#45).

usage: reads_spike_plan.py COCCI_TRUTH.tsv COCCI_R1.list TITRATION_DIR OUT_DIR [N_COCCI=8] [N_ZYGO=12] [FRACS=1,5,10,20,30,50]
Hosts: N_COCCI Coccidioides read sets with a species truth (half C. immitis, half C. posadasii; species is in the reference, so the
species call can be followed) and N_ZYGO ZyGoLife read sets whose genus is in the reference (genus truth from the label, as in the titration).
Writes OUT_DIR/hosts.tsv (id, r1, r2) and OUT_DIR/truth.tsv (id, phylum..species, host, donor, frac for every variant, including the host alone).
"""
import csv, os, random, sys

cocci_truth, cocci_list, titr, out = sys.argv[1:5]
nc = int(sys.argv[5]) if len(sys.argv) > 5 else 8
nz = int(sys.argv[6]) if len(sys.argv) > 6 else 12
fracs = [int(x) for x in (sys.argv[7] if len(sys.argv) > 7 else "1,5,10,20,30,50").split(",")]
rng = random.Random(11)
os.makedirs(out, exist_ok=True)
R = ["phylum", "class", "order", "family", "genus", "species"]
paths = {os.path.basename(p.strip()).rsplit(".fastq", 1)[0]: p.strip() for p in open(cocci_list) if p.strip()}
cands = {"Coccidioides immitis": [], "Coccidioides posadasii": []}
for r in csv.DictReader(open(cocci_truth), delimiter="\t"):
    p = paths.get(r["id"])
    if not p or r["species"] not in cands or "_R1_" not in p: continue
    p2 = p.replace("_R1_", "_R2_")
    if os.path.exists(p) and os.path.exists(p2) and os.path.getsize(p) > 60_000_000: cands[r["species"]].append((r, p, p2))
hosts = []
for sp, lst in cands.items():
    rng.shuffle(lst)
    for r, p, p2 in lst[: nc // 2]:
        hosts.append((r["id"], p, p2, [r[k] for k in R]))
sheet = {}
for r in csv.DictReader(open(f"{titr}/sheet.tsv"), delimiter="\t"):
    sheet.setdefault(r["id"], []).append(r["path"])
z = [r for r in csv.DictReader(open(f"{titr}/truth.tsv"), delimiter="\t") if r.get("genus_in_ref") == "True" and len(sheet.get(r["id"], [])) == 2]
rng.shuffle(z)
for r in z[:nz]:
    a, b = sheet[r["id"]]
    hosts.append((r["id"], a, b, [r[k] for k in R]))
with open(f"{out}/hosts.tsv", "w") as fh:
    for i, a, b, _ in hosts: fh.write(f"{i}\t{a}\t{b}\n")
with open(f"{out}/truth.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["id"] + R + ["host", "donor", "frac"])
    for i, _, _, lin in hosts:
        w.writerow([f"{i}__none_0"] + lin + [i, "none", 0])
        for d in ("bacteria", "human"):
            for f in fracs: w.writerow([f"{i}__{d}_{f}"] + lin + [i, d, f])
print(len(hosts), "hosts", file=sys.stderr)
