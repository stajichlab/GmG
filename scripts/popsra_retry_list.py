#!/usr/bin/env python3
"""Which runs of a manifest need another go, from the result lines `popseq fetch` wrote.

usage: popsra_retry_list.py MANIFEST RESULTS [-o RETRY.csv]

MANIFEST is the CSV or TSV the chunk jobs read (a run_accession column); RESULTS is the tab-separated results file the fetches
appended to (header: run tier route layout ... verified outcome [sampling fallback_from]). The LAST line of a run decides its
state, and the policy follows PopSeq spec 002's exit codes:

  ok                   done
  transient            retry (the transfer or a lookup failed: a retry is reasonable)
  verify_failed        retry once; a run with two verify_failed lines is set aside ("quarantine"), as a repeat means the
                       source, not the transfer, is bad
  hook_failed          set aside ("inspect"): the download is ok, the sketch check or the sketcher failed; a person looks
  unavailable_at_ncbi  final, not retried (no archive serves the run)
  unpaired             final, not retried (pairs were required and no archive serves them)
  (no line)            retry: the task was killed or never ran

The retry manifest keeps the manifest's rows and columns, in manifest order. A summary of every state is printed to stderr.
"""
import argparse
import csv
import sys
from collections import Counter

FINAL = {"unavailable_at_ncbi", "unpaired"}


def read_table(path):
    with open(path, newline="") as fh:
        head = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in head else ","
        return list(csv.DictReader(fh, delimiter=delim)), delim


def classify(runs, results):
    """`{run: state}` with state in done, retry, quarantine, inspect, final."""
    lines = {}
    for r in results:
        lines.setdefault(r["run"], []).append(r["outcome"])
    out = {}
    for run in runs:
        seen = lines.get(run, [])
        last = seen[-1] if seen else ""
        if last == "ok":
            out[run] = "done"
        elif last in FINAL:
            out[run] = "final"
        elif last == "hook_failed":
            out[run] = "inspect"
        elif last == "verify_failed":
            out[run] = "quarantine" if seen.count("verify_failed") >= 2 else "retry"
        else:  # transient, or no line at all, or an outcome this script does not know
            out[run] = "retry"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest")
    ap.add_argument("results")
    ap.add_argument("-o", "--out", help="write the retry manifest here (default: stdout)")
    a = ap.parse_args(argv)
    rows, delim = read_table(a.manifest)
    if not rows or "run_accession" not in rows[0]:
        print(f"{a.manifest}: no run_accession column", file=sys.stderr)
        return 2
    results, _ = read_table(a.results)
    runs = [r["run_accession"] for r in rows]
    state = classify(runs, results)
    summary = Counter(state.values())
    print("state\truns", file=sys.stderr)
    for k in ("done", "retry", "quarantine", "inspect", "final"):
        print(f"{k}\t{summary.get(k, 0)}", file=sys.stderr)
    todo = [r for r in rows if state[r["run_accession"]] == "retry"]
    out = open(a.out, "w", newline="") if a.out else sys.stdout
    w = csv.DictWriter(out, fieldnames=list(rows[0].keys()), delimiter=delim, lineterminator="\n")
    w.writeheader()
    w.writerows(todo)
    if a.out:
        out.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
