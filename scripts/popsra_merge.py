#!/usr/bin/env python3
"""Join the per-task TSV files of a chunk-job run into one file each, in $WORK/merged/.

usage: popsra_merge.py WORK

Each array task of popsra_fetch_chunks.slurm writes results.tsv, profile.tsv, exits.tsv, hook_times.tsv and sampling_check.tsv to
WORK/tasks/<i>/ (appends from several nodes to one shared file can interleave or lose lines). The inputs for each name are
WORK/<name>.tsv if it exists (the flat layout of the earlier runs) followed by the task files in numeric task order. The first
line of the first file is kept as the header and a later line equal to it is dropped (a task that was run twice). The output goes
to a temporary name and is renamed; nothing is deleted, and the inputs are never modified, so the merge can be repeated.
"""
import os
import re
import sys

NAMES = ("results", "profile", "exits", "hook_times", "sampling_check")
# exits.tsv and hook_times.tsv have no header line; the others do (their first line starts with "run" or "id")
HEADERED = {"results": "run", "profile": "run", "sampling_check": "id"}


def task_files(work, name):
    tasks = os.path.join(work, "tasks")
    if not os.path.isdir(tasks):
        return []
    found = []
    for d in os.listdir(tasks):
        f = os.path.join(tasks, d, f"{name}.tsv")
        if re.fullmatch(r"\d+", d) and os.path.isfile(f):
            found.append((int(d), f))
    return [f for _, f in sorted(found)]


def merge(work, name):
    inputs = []
    flat = os.path.join(work, f"{name}.tsv")
    if os.path.isfile(flat):
        inputs.append(flat)
    inputs += task_files(work, name)
    if not inputs:
        return None
    first_col = HEADERED.get(name)
    header = None
    lines = []
    for path in inputs:
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip():
                    continue  # a blank line is never a record
                is_header = first_col is not None and line.split("\t")[0] == first_col
                if is_header:
                    if header is None:
                        header = line
                    continue  # the header of a later file, or of a task run twice
                lines.append(line)
    out_dir = os.path.join(work, "merged")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{name}.tsv")
    tmp = f"{out}.tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        if header is not None:
            fh.write(header + "\n")
        for line in lines:
            fh.write(line + "\n")
    os.replace(tmp, out)
    return out, len(inputs), len(lines)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1 or not os.path.isdir(argv[0]):
        print(__doc__, file=sys.stderr)
        return 2
    for name in NAMES:
        r = merge(argv[0], name)
        if r:
            print(f"{r[0]}: {r[2]} lines from {r[1]} files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
