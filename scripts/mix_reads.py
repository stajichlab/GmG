#!/usr/bin/env python3
"""Mix real read pairs of two or three samples (#45): the first N * pct / 100 pairs of each parent, written as plain paired FASTQ.

usage: mix_reads.py OUT_DIR NAME N PCTS A_R1 A_R2 [B_R1 B_R2 [C_R1 C_R2]]
  PCTS  comma-separated shares of the N pairs, one per parent, summing to 100 (e.g. 90,10 or 60,30,10)
Writes OUT_DIR/NAME_R1.fastq and OUT_DIR/NAME_R2.fastq and prints the two sample-sheet rows (id, path). A share of 0 skips that parent
(a pure control is `100,0`). Each parent is read from its start, so the parents' reads are different sequences. Plain or .gz input;
read-only on the inputs. Fails if a parent has fewer pairs than its share needs.
"""
import gzip, itertools, os, sys

out, name, n, pcts = sys.argv[1:5]
n = int(n)
pcts = [float(x) for x in pcts.split(",")]
files = sys.argv[5:]
assert len(files) == 2 * len(pcts), "one R1 and one R2 per share"
assert abs(sum(pcts) - 100) < 1e-6, "shares must sum to 100"
os.makedirs(out, exist_ok=True)
counts = [round(n * p / 100) for p in pcts]
counts[0] += n - sum(counts)   # rounding goes to the first (largest) parent


def reads(path, pairs):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        got = list(itertools.islice(fh, pairs * 4))
    if len(got) < pairs * 4:
        sys.exit(f"{path} has {len(got) // 4} read pairs, {pairs} needed")
    return got


for k, tag in (("R1", 0), ("R2", 1)):
    with open(os.path.join(out, f"{name}_{k}.fastq"), "w") as w:
        for i, c in enumerate(counts):
            if c:
                w.writelines(reads(files[2 * i + tag], c))
for k in ("R1", "R2"):
    print(f"{name}\t{os.path.join(out, f'{name}_{k}.fastq')}")
