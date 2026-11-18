#!/usr/bin/env python3
"""The species name each run of a manifest is filed under, as a table for the species check in popsra_sketch_hook.sh.

usage: popsra_filed_table.py MANIFEST -o OUT.tsv

MANIFEST is a CSV or TSV with a run column (`run_accession`, or `run`) and a name column: `effective_species` (curated, else label,
else NCBI: PopSeq spec 009) when the manifest has it, otherwise `species_name`. OUT is `id<TAB>filed` with a header line, written to a
temporary name and renamed; a row with no name is left out (the hook records such a run as unchecked). Exit 2 for an unusable manifest.
"""
import argparse
import csv
import os
import sys


def filed_names(path):
    with open(path, newline="") as fh:
        head = fh.readline()
        fh.seek(0)
        delim = "\t" if head.count("\t") > head.count(",") else ","
        rd = csv.DictReader(fh, delimiter=delim)
        cols = rd.fieldnames or []
        run_col = next((c for c in ("run_accession", "run") if c in cols), None)
        name_col = next((c for c in ("effective_species", "species_name") if c in cols), None)
        if run_col is None or name_col is None:
            raise ValueError(f"{path}: needs a run column (run_accession or run) and a name column (effective_species or species_name); has {cols}")
        out, seen = [], set()
        for r in rd:
            run, name = (r[run_col] or "").strip(), (r[name_col] or "").strip()
            if run and name and run not in seen:
                seen.add(run)
                out.append((run, name.replace("\t", " ")))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("manifest")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args(argv)
    try:
        rows = filed_names(a.manifest)
    except (OSError, ValueError) as e:
        print(f"popsra_filed_table: {e}", file=sys.stderr)
        return 2
    tmp = f"{a.out}.tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write("id\tfiled\n")
        for run, name in rows:
            fh.write(f"{run}\t{name}\n")
    os.replace(tmp, a.out)
    print(f"{len(rows)} runs with a filed name -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
