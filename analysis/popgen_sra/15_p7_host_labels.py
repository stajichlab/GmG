#!/usr/bin/env python3
"""#155: labels for the P7 P. oryzae set from the PopSeq catalog (no new query: catalog/runs/fungi_all_runs.csv.gz).
usage: 15_p7_host_labels.py CATALOG_DIR MANIFEST SKETCH_DIR OUT.tsv
Writes id, host, study, country for every sketched run in the manifest. host = the catalog's host_name when its host_status is ok, else blank;
rare hosts (fewer than MIN_GROUP runs) are pooled into 'other_host' so the label test has groups it can use."""
import csv, gzip, collections, os, sys
cat, manifest, sk, out = sys.argv[1:5]
MIN_GROUP = 5
ids = {f[:-5] for f in os.listdir(sk) if f.endswith(".json") and ".half" not in f}
rows = {r["run_accession"]: r for r in csv.DictReader(gzip.open(f"{cat}/runs/fungi_all_runs.csv.gz", "rt"))}
man = {r["run_accession"] for r in csv.DictReader(open(manifest))}
host = {i: (rows[i]["host_name"] if rows[i].get("host_status") == "ok" else "") for i in ids if i in rows and i in man}
n = collections.Counter(host.values())
with open(out, "w") as f:
    f.write("id\thost\tstudy\tcountry\n")
    for i in sorted(host):
        h = host[i] if host[i] and n[host[i]] >= MIN_GROUP else ("other_host" if host[i] else "")
        f.write(f"{i}\t{h}\t{rows[i]['study_accession']}\t{rows[i]['country_name']}\n")
print(len(host), "labelled;", dict(collections.Counter(("other_host" if v and n[v] < MIN_GROUP else v) for v in host.values())))
