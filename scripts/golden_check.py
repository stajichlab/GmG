#!/usr/bin/env python3
"""Golden-panel regression test (#49): real genomes with pinned inputs and pinned expected calls.

usage: golden_check.py --gmg BIN --index INDEX.gmgi [--cache DIR] [--update] [--expected FILE] [--threads N]
  data/golden/panel.tsv       the 17 genomes: accession, NCBI FTP URL and SHA-256 of the file (an assembly version never changes)
  data/golden/expected.tsv    what `gmg run` must give with the pinned index (data/golden/INDEX.txt names it)
Downloads missing files into --cache (default ~/.cache/gmg-golden, checked against the SHA-256), runs `gmg run`, and compares per genome:
status, best_rank, the call at every rank, the QC flags, the top hit and the ANI (to 0.0005). Exit code 1 and a diff when anything differs.
--update rewrites expected.tsv from the current run (do this only for an intended change, and say why in the commit message).
"""
import argparse, csv, hashlib, os, subprocess, sys, tempfile, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
GOLD = os.path.join(HERE, "..", "data", "golden")
COLS = ["status", "best_rank", "call_phylum", "call_class", "call_order", "call_family", "call_genus", "call_species", "qc_flags", "top_hit"]
ANI_TOL = 0.0005

ap = argparse.ArgumentParser()
ap.add_argument("--gmg", default="gmg")
ap.add_argument("--index", default="")
ap.add_argument("--cache", default=os.path.expanduser("~/.cache/gmg-golden"))
ap.add_argument("--expected", default=os.path.join(GOLD, "expected.tsv"))
ap.add_argument("--update", action="store_true")
ap.add_argument("--threads", default="4")
ap.add_argument("--fetch-only", action="store_true", help="download and verify the panel into --cache, then stop")
a = ap.parse_args()

panel = list(csv.DictReader(open(os.path.join(GOLD, "panel.tsv")), delimiter="\t"))
os.makedirs(a.cache, exist_ok=True)
work = tempfile.mkdtemp(prefix="gmg_golden_")
os.makedirs(f"{work}/in")
for p in panel:
    f = f"{a.cache}/{p['id']}.fna.gz"
    if not (os.path.exists(f) and hashlib.sha256(open(f, "rb").read()).hexdigest() == p["sha256"]):
        print("fetching", p["id"], file=sys.stderr)
        data = urllib.request.urlopen(p["url"], timeout=600).read()
        if hashlib.sha256(data).hexdigest() != p["sha256"]:
            sys.exit(f"checksum mismatch for {p['id']}: the file at NCBI changed")
        open(f, "wb").write(data)
    os.symlink(f, f"{work}/in/{p['id']}.fna.gz")
if a.fetch_only:
    print(f"{len(panel)} genomes verified in {a.cache}")
    sys.exit(0)
r = subprocess.run([a.gmg, "--threads", a.threads, "-q", "run", f"{work}/in", "--index", a.index, "-o", f"{work}/out", "--no-report"],
                   capture_output=True, text=True)
if r.returncode != 0:
    sys.exit("gmg run failed:\n" + r.stderr)
got = {x["query_id"]: x for x in csv.DictReader(open(f"{work}/out/results.tsv"), delimiter="\t")}


def norm(x):
    d = {c: x[c] for c in COLS}
    d["qc_flags"] = ",".join(sorted(f for f in x["qc_flags"].split(",") if f))
    d["ani"] = x["ani"]
    return d


if a.update:
    with open(a.expected, "w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["id"] + COLS + ["ani"])
        for p in panel:
            d = norm(got[p["id"]])
            w.writerow([p["id"]] + [d[c] for c in COLS] + [d["ani"]])
    print("expected.tsv rewritten from this run")
    sys.exit(0)

exp = {x["id"]: x for x in csv.DictReader(open(a.expected), delimiter="\t")}
bad = []
for p in panel:
    i = p["id"]
    if i not in got:
        bad.append(f"{i}: no result"); continue
    g, e = norm(got[i]), exp[i]
    for c in COLS:
        if g[c] != e[c]:
            bad.append(f"{i} {c}: expected {e[c]!r}, got {g[c]!r}")
    ga, ea = g["ani"], e["ani"]
    if (ga == "") != (ea == "") or (ga != "" and abs(float(ga) - float(ea)) > ANI_TOL):
        bad.append(f"{i} ani: expected {ea}, got {ga}")
if bad:
    print(f"golden panel: {len(bad)} difference(s) in {len({b.split(':')[0].split(' ')[0] for b in bad})} of {len(panel)} genomes", file=sys.stderr)
    print("\n".join(bad), file=sys.stderr)
    sys.exit(1)
print(f"golden panel: {len(panel)} genomes identical to {os.path.relpath(a.expected)}")
