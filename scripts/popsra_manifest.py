#!/usr/bin/env python3
"""Print the run_accession of every data row of a manifest, one per line (blank rows skipped).

usage: popsra_manifest.py MANIFEST

The manifest is a CSV or TSV with a run_accession column, as `popseq fetch` reads it. It is read with Python's csv module, so a
quoted field that holds a comma (a `geo_loc_name` such as "USA: California, Riverside") does not shift the columns. Exit 2 if
the file cannot be read or has no run_accession column; a manifest with a header and no rows prints nothing and exits 0.
"""
import csv
import sys


def runs(path):
    with open(path, newline="") as fh:
        head = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in head else ","
        reader = csv.DictReader(fh, delimiter=delim)
        if not reader.fieldnames or "run_accession" not in reader.fieldnames:
            raise ValueError(f"no run_accession column in {path}")
        return [r["run_accession"].strip() for r in reader if (r.get("run_accession") or "").strip()]


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        print(__doc__, file=sys.stderr)
        return 2
    try:
        found = runs(argv[0])
    except (OSError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 2
    for r in found:
        print(r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
