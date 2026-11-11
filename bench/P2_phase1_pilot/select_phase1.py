#!/usr/bin/env python3
"""P2: choose N isolates of one species from the PopSeq panel, spread across studies, for the phase-1 pilot.

usage: select_phase1.py PANEL_METADATA.csv TAXID N OUT_MANIFEST.csv [MIN_DEPTH MAX_DEPTH]
Round robin over studies (sorted by name): the first remaining run of each study in turn, until N are
chosen, so the set is deterministic and not one study's batch. Only runs with an estimated depth between
MIN_DEPTH and MAX_DEPTH (default 25 and 200) are considered: the coverage tier keeps 25x but the transfer is the
whole run, so a 600x run costs a 24-fold larger download than it keeps. Writes a manifest with a run_accession column."""
import csv
import sys
from itertools import zip_longest

meta, taxid, n, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
lo, hi = (float(sys.argv[5]), float(sys.argv[6])) if len(sys.argv) > 6 else (25.0, 200.0)
by_study = {}
for r in csv.DictReader(open(meta)):
    if r["species_taxid"] == taxid and r["est_depth_x"] and lo <= float(r["est_depth_x"]) <= hi:
        by_study.setdefault(r["study_accession"], []).append(r)
for rows in by_study.values():
    rows.sort(key=lambda r: r["run_accession"])
order = [by_study[s] for s in sorted(by_study)]
chosen = [r for group in zip_longest(*order) for r in group if r][:n]
fields = ["run_accession", "species_name", "study_accession", "country", "est_depth_x", "base_count"]
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
    w.writeheader()
    w.writerows(chosen)
print(f"{len(chosen)} runs from {len({r['study_accession'] for r in chosen})} studies of {len(by_study)} for taxid {taxid}")
tot = sum(int(r["base_count"]) for r in chosen)
for r in chosen:
    print(" ", r["run_accession"], r["study_accession"], r["est_depth_x"] + "x", int(r["base_count"]) // 10**6, "Mb")
print(f"total {tot / 1e6:.0f} Mb of bases; the transfer is about {tot * 0.28 / 1e9:.1f} to {tot * 0.65 / 1e9:.1f} GB "
      "(0.28 to 0.65 bytes per base, the `.sra` range measured so far, n = 3)")
