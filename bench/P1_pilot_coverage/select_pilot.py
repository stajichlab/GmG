#!/usr/bin/env python3
"""P1: choose a few runs of a few species from the PopSeq panel and write a plan CSV for `popseq fetch`.

usage: select_pilot.py PANEL_METADATA.csv GENOME_SIZES.csv OUT_PLAN.csv
Two isolates for each of three small-genome species (taxids below), each a run of 35-70x and
300-650 Mb of bases (so a 25x target keeps about half of it), the smallest such runs first.
The plan carries only ENA values (md5, bytes, URLs) and no NCBI columns: with no NCBI counts the
SRA routes fail at once ("cannot be verified") and ENA, the only route a laptop without the SRA
toolkit can use, does the work. Needs internet for the ENA lookup.
"""
import csv
import sys

sys.path.insert(0, "popseq")
from popseq.lookup import ena_lookup  # noqa: E402

SPECIES = ["5478", "5007", "5476"]   # N. glabratus, B. bruxellensis, C. albicans
PER_SPECIES = 2

meta_path, gs_path, out = sys.argv[1:4]
gs = {r["species_taxid"]: float(r["genome_mb"]) for r in csv.DictReader(open(gs_path)) if r["genome_mb"]}
chosen = []
for tid in SPECIES:
    rows = [r for r in csv.DictReader(open(meta_path))
            if r["species_taxid"] == tid and r["est_depth_x"] and 35 <= float(r["est_depth_x"]) <= 70
            and r["base_count"] and 300e6 <= int(r["base_count"]) <= 650e6]
    rows.sort(key=lambda r: int(r["base_count"]))
    chosen += rows[:PER_SPECIES]
ena = ena_lookup([r["run_accession"] for r in chosen])
fields = ["run", "taxid", "species", "study", "read_count", "base_count", "fastq_bytes", "fastq_md5",
          "fastq_ftp", "fastq_aspera", "role", "ncbi_spots", "ncbi_spots_with_mates"]
total = 0
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    for r in chosen:
        e = ena[r["run_accession"]]
        b = sum(int(x) for x in e["fastq_bytes"].split(";") if x)
        total += b
        w.writerow({"run": r["run_accession"], "taxid": r["species_taxid"], "species": r["species_name"],
                    "study": r["study_accession"], "read_count": e["read_count"], "base_count": e["base_count"],
                    "fastq_bytes": e["fastq_bytes"], "fastq_md5": e["fastq_md5"], "fastq_ftp": e["fastq_ftp"],
                    "fastq_aspera": e["fastq_aspera"], "role": "pilot"})
        print(f"{r['run_accession']}  {r['species_name']:28s} est {float(r['est_depth_x']):5.1f}x  "
              f"{int(r['base_count'])/1e6:5.0f} Mb  fastq {b/1e6:6.1f} MB  {gs[r['species_taxid']]:.2f} Mb genome")
print(f"{len(chosen)} runs, {total/1e6:.0f} MB of FASTQ to download; at 0.2 MB/s that is about {total/2e5/60:.0f} min")
