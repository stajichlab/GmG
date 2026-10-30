#!/usr/bin/env python3
"""Add a deterministic `split` column (cal | eval) to a truth table (#32).

usage: reads_calibration_split.py TRUTH.tsv OUT.tsv [--id-column id]
The split depends only on the sample id: the first byte of SHA-1(id) is even -> `cal` (the calibration set: the cutoff was chosen on these),
odd -> `eval` (the evaluation set: the numbers that are reported). The same id always lands in the same half, whatever else is in the table,
so a result table can be re-split without re-running anything.
"""
import argparse, csv, hashlib

ap = argparse.ArgumentParser()
ap.add_argument("truth"); ap.add_argument("out"); ap.add_argument("--id-column", default="id")
a = ap.parse_args()


def split_of(sample_id):
    return "cal" if hashlib.sha1(sample_id.encode()).digest()[0] % 2 == 0 else "eval"


rows = list(csv.DictReader(open(a.truth), delimiter="\t"))
cols = list(rows[0].keys()) + (["split"] if "split" not in rows[0] else [])
with open(a.out, "w", newline="") as o:
    w = csv.DictWriter(o, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        r["split"] = split_of(r[a.id_column])
        w.writerow(r)
n_cal = sum(r["split"] == "cal" for r in rows)
print(f"{len(rows)} samples: {n_cal} calibration, {len(rows) - n_cal} evaluation")
