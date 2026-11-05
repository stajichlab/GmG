#!/usr/bin/env python3
"""02_ena_urls.py : add the ENA FASTQ addresses and checksums to each species' runs.tsv (-> runs_ftp.tsv).
usage: 02_ena_urls.py DIR   (DIR/<slug>/runs.tsv; needs network)  Standard library only."""
import csv, glob, os, sys, urllib.parse, urllib.request, time
API = "https://www.ebi.ac.uk/ena/portal/api/search"
for path in sorted(glob.glob(os.path.join(sys.argv[1], "*", "runs.tsv"))):
    rows = list(csv.DictReader(open(path, newline=""), delimiter="\t"))
    info = {}
    for i in range(0, len(rows), 50):
        batch = rows[i:i + 50]
        q = " OR ".join(f"run_accession={r['run_accession']}" for r in batch)
        data = urllib.parse.urlencode({"result": "read_run", "query": q, "fields": "run_accession,fastq_ftp,fastq_md5", "format": "tsv", "limit": 0}).encode()
        for attempt in range(5):
            try:
                txt = urllib.request.urlopen(urllib.request.Request(API, data=data), timeout=120).read().decode()
                break
            except Exception as e:
                if attempt == 4: raise
                time.sleep(5 * (attempt + 1))
        for line in txt.strip().split("\n")[1:]:
            p = line.split("\t")
            if len(p) >= 3: info[p[0]] = (p[1], p[2])
    out = os.path.join(os.path.dirname(path), "runs_ftp.tsv")
    n = 0
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(list(rows[0].keys()) + ["fastq_ftp", "fastq_md5"])
        for r in rows:
            ftp, md5 = info.get(r["run_accession"], ("", ""))
            if len(ftp.split(";")) < 2: continue  # needs both mates
            w.writerow(list(r.values()) + [ftp, md5]); n += 1
    print(f"{os.path.basename(os.path.dirname(path))}: {n} of {len(rows)} runs with paired FASTQ on ENA")
