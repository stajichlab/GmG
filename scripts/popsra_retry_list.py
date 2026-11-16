#!/usr/bin/env python3
"""Which runs of a manifest need another go, from the result lines `popseq fetch` wrote.

usage: popsra_retry_list.py MANIFEST RESULTS [RESULTS ...] [-o RETRY.csv] [--max-attempts N] [--exit-status]

MANIFEST is the CSV or TSV the chunk jobs read (a run_accession column); each RESULTS is a tab-separated results file the
fetches appended to (header: run tier route layout ... verified outcome [sampling fallback_from]); give the merged file
(`popsra_merge.py`) or the per-task files (`$WORK/tasks/*/results.tsv`). A missing RESULTS file counts as no lines. The LAST line
of a run decides its state, and the policy follows PopSeq spec 002's exit codes:

  ok                   done
  transient            retry, up to --max-attempts lines for the run (default 3); then set aside ("inspect"): a transfer
                       that fails three times is a person's problem, not a loop
  verify_failed        retry once; a run with two verify_failed lines is set aside ("quarantine"), as a repeat means the
                       source, not the transfer, is bad
  hook_failed          set aside ("inspect"): the download is ok, the sketch check or the sketcher failed; a person looks
  unavailable_at_ncbi  final, not retried (no archive serves the run)
  unpaired             final, not retried (pairs were required and no archive serves them)
  an outcome not in this list   set aside ("inspect"): a new outcome from popseq must be read by a person, never retried
                       without limit
  (no line)            retry: the task was killed or never ran

The retry manifest keeps the manifest's rows and columns, in manifest order; a manifest with no rows left gives a header-only
file, which the chunk jobs read as zero chunks. A summary of every state is printed to stderr. Exit 0, unless --exit-status:
then 10 when runs remain to retry and 0 when none do, so a driver can loop until done. Exit 2 for a manifest with no
run_accession column.
"""
import argparse
import csv
import sys
from collections import Counter

FINAL = {"unavailable_at_ncbi", "unpaired"}
KNOWN = {"ok", "transient", "verify_failed", "hook_failed"} | FINAL
DEFAULT_MAX_ATTEMPTS = 3


def read_table(path):
    with open(path, newline="") as fh:
        head = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in head else ","
        return list(csv.DictReader(fh, delimiter=delim)), delim


# the columns of a result line (popseq.fetch.FETCH_FIELDS): used for a results file whose header line is missing, as happens when
# two fetches race to create the file
DEFAULT_FIELDS = [
    "run", "tier", "route", "layout", "bases_kept", "reads_kept", "coverage_target", "coverage_reached",
    "seconds", "verified", "outcome", "sampling", "fallback_from",
]


def read_results(paths):
    """The result rows of every existing file; a missing file is no lines. A file whose first line is not a header is read with
    the default columns (no line is lost); a repeated header, a blank or a short line is skipped and counted."""
    rows, skipped = [], 0
    for path in paths:
        try:
            with open(path, newline="") as fh:
                text = [l for l in fh.read().splitlines() if l.strip()]
        except FileNotFoundError:
            continue
        if not text:
            continue
        delim = "\t" if "\t" in text[0] else ","
        fields = None if text[0].split(delim)[0] == "run" else DEFAULT_FIELDS
        for rec in csv.reader(text, delimiter=delim):
            if fields is None:
                fields = rec
                continue
            if rec and rec[0] == "run":  # a header repeated in the middle
                skipped += 1
                continue
            r = dict(zip(fields, rec))
            if r.get("run") and r.get("outcome"):
                rows.append(r)
            else:
                skipped += 1
    return rows, skipped


def classify(runs, results, max_attempts=DEFAULT_MAX_ATTEMPTS):
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
        elif last == "transient":
            out[run] = "inspect" if len(seen) >= max_attempts else "retry"
        elif last == "":  # no line: the task was killed or never ran
            out[run] = "retry"
        else:  # an outcome this script does not know
            out[run] = "inspect"
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest")
    ap.add_argument("results", nargs="+")
    ap.add_argument("-o", "--out", help="write the retry manifest here (default: stdout)")
    ap.add_argument("--max-attempts", type=int, default=DEFAULT_MAX_ATTEMPTS, help="result lines after which a transient run is set aside (default 3)")
    ap.add_argument("--exit-status", action="store_true", help="exit 10 when runs remain to retry, 0 when none do")
    a = ap.parse_args(argv)
    try:
        rows, delim = read_table(a.manifest)
        with open(a.manifest, newline="") as fh:
            fields = next(csv.reader(fh, delimiter=delim), [])
    except OSError as exc:
        print(f"{a.manifest}: {exc}", file=sys.stderr)
        return 2
    if "run_accession" not in fields:
        print(f"{a.manifest}: no run_accession column", file=sys.stderr)
        return 2
    results, skipped = read_results(a.results)
    if skipped:
        print(f"note: {skipped} line(s) of the results were not result rows and were skipped", file=sys.stderr)
    runs = [r["run_accession"] for r in rows]
    state = classify(runs, results, a.max_attempts)
    summary = Counter(state.values())
    print("state\truns", file=sys.stderr)
    for k in ("done", "retry", "quarantine", "inspect", "final"):
        print(f"{k}\t{summary.get(k, 0)}", file=sys.stderr)
    todo = [r for r in rows if state[r["run_accession"]] == "retry"]
    out = open(a.out, "w", newline="") if a.out else sys.stdout
    w = csv.DictWriter(out, fieldnames=fields, delimiter=delim, lineterminator="\n", extrasaction="ignore")
    w.writeheader()
    w.writerows(todo)
    if a.out:
        out.close()
    return 10 if (a.exit_status and todo) else 0


if __name__ == "__main__":
    raise SystemExit(main())
