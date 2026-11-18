#!/usr/bin/env python3
"""Compare two species-match runs of the same sketches made with different indexes (for example bfd_s1000 against bfd_s10000, the index the fetch hook would load).
usage: compare_indexes.py DIR_A DIR_B   (each holds sm_<set>.tsv from `popsketch species-match`)
Prints, over all runs present in both: the status agreement, the agreement of exclude_default, and every run whose status differs, with both statuses and purities.
A run whose status differs but whose exclude_default is the same is a difference in wording only; a difference in exclude_default changes a panel."""
import csv, glob, os, sys, collections
a_dir, b_dir = sys.argv[1:3]
def load(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "sm_*.tsv"))):
        for r in csv.DictReader(open(f), delimiter="\t"):
            out[r["id"]] = r
    return out
A, B = load(a_dir), load(b_dir)
ids = sorted(set(A) & set(B))
print(f"{len(A)} runs in A, {len(B)} in B, {len(ids)} in both")
same = sum(A[i]["status"] == B[i]["status"] for i in ids)
same_x = sum(A[i]["exclude_default"] == B[i]["exclude_default"] for i in ids)
print(f"status identical: {same}/{len(ids)} ({100*same/len(ids):.1f}%); exclude_default identical: {same_x}/{len(ids)}")
conf = collections.Counter((A[i]["status"], B[i]["status"]) for i in ids if A[i]["status"] != B[i]["status"])
print("status changes A -> B:")
for (x, y), n in sorted(conf.items(), key=lambda kv: -kv[1]):
    print(f"  {x:24s} -> {y:24s} {n}")
print("runs whose exclude_default differs:")
for i in ids:
    if A[i]["exclude_default"] != B[i]["exclude_default"]:
        print(f"  {i}  {A[i]['filed_name']}  A={A[i]['status']}({A[i]['purity']})  B={B[i]['status']}({B[i]['purity']})")
