#!/usr/bin/env python3
"""#159 step 4: estimated genome size against the species' assembly size, per run, from finished popsketch sets.
usage (from popsra): 12_size_vs_assembly.py OUTSUFFIX SLUG [SLUG ...]   (qc from <slug>_<OUTSUFFIX>/qc.tsv, assembly size = genome_mb in <slug>/runs.tsv, PopSeq catalog)
Prints one TSV (slug, run, species, assembly_mb, est_mb, ratio, coverage_basis, status) to stdout and a per-species summary to stderr.
A run is an outlier when est/assembly is outside 0.75-1.33 (the second-genome and contaminant cases land far outside; real size variation within a species does not)."""
import csv, statistics, sys
suffix, slugs = sys.argv[1], sys.argv[2:]
print("slug\trun\tspecies\tassembly_mb\test_mb\tratio\tcoverage_basis\tstatus\toutlier")
for s in slugs:
    runs = {r["run_accession"]: r for r in csv.DictReader(open(f"{s}/runs.tsv"), delimiter="\t")}
    n = out = 0; ratios = []
    for q in csv.DictReader(open(f"{s}_{suffix}/qc.tsv"), delimiter="\t"):
        r = runs.get(q["id"])
        if not r or not q["est_size_mb"] or not r.get("genome_mb"): continue
        a, e = float(r["genome_mb"]), float(q["est_size_mb"])
        ratio = e / a; o = not 0.75 <= ratio <= 1.33
        n += 1; out += o; ratios.append(ratio)
        print(f'{s}\t{q["id"]}\t{r["species"]}\t{a}\t{e}\t{ratio:.2f}\t{q["coverage_basis"]}\t{q["status"]}\t{"outlier" if o else ""}')
    print(f"{s}: {n} runs with a size, median ratio {statistics.median(ratios):.2f}, {out} outside 0.75-1.33 ({100*out/n:.1f}%)", file=sys.stderr)
