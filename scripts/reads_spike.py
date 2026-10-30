#!/usr/bin/env python3
"""Contaminated read sets (#45): real fungal read pairs plus simulated bacterial or human read pairs, written as plain FASTQ.

usage: reads_spike.py ID HOST_R1 HOST_R2 BACTERIA.fna HUMAN.fna OUT_DIR [N_PAIRS=300000] [FRACS=1,5,10,20,30,50] [SEED=3]
Takes the first N_PAIRS read pairs of the host and writes, in OUT_DIR, the host alone (`ID__none_0_R1.fastq`, `_R2.fastq`) and, for each donor
(bacteria, human) and each share f in FRACS (percent of all read pairs), the host pairs plus M = N * f / (100 - f) simulated pairs from the donor genome:
2 x 150 bp, fragment 300 bp, error-free, uniform over the genome (both strands), quality `F`. Names: `ID__bacteria_10_R1.fastq`. Prints the
sample sheet rows (id, path) of every variant to stdout. Read-only on the host files.
"""
import gzip, random, sys, os

sid, r1, r2, bac, hum, out = sys.argv[1:7]
n = int(sys.argv[7]) if len(sys.argv) > 7 else 300_000
fracs = [int(x) for x in (sys.argv[8] if len(sys.argv) > 8 else "1,5,10,20,30,50").split(",")]
rng = random.Random(int(sys.argv[9]) if len(sys.argv) > 9 else 3)
os.makedirs(out, exist_ok=True)
COMP = str.maketrans("ACGT", "TGCA")


def head(path, pairs):
    op = gzip.open if path.endswith(".gz") else open
    lines = []
    with op(path, "rt") as fh:
        for i, l in enumerate(fh):
            if i >= pairs * 4: break
            lines.append(l)
    return "".join(lines)


def genome(path):
    seq = []
    for l in open(path):
        if not l.startswith(">"): seq.append(l.strip().upper())
    return "".join(seq)


def sim(g, m, tag):
    a, b = [], []
    L = len(g)
    q = "F" * 150
    for i in range(m):
        while True:
            s = rng.randrange(0, L - 300)
            f = g[s:s + 300]
            if "N" not in f: break
        if rng.random() < 0.5: f = f.translate(COMP)[::-1]
        a.append(f"@sim_{tag}_{i}/1\n{f[:150]}\n+\n{q}\n")
        b.append(f"@sim_{tag}_{i}/2\n{f[150:].translate(COMP)[::-1]}\n+\n{q}\n")
    return "".join(a), "".join(b)


h1, h2 = head(r1, n), head(r2, n)
rows = []


def emit(name, x1, x2):
    for k, x in (("R1", x1), ("R2", x2)):
        p = os.path.join(out, f"{name}_{k}.fastq")
        open(p, "w").write(x)
        rows.append(f"{name}\t{p}")


emit(f"{sid}__none_0", h1, h2)
for donor, path in (("bacteria", bac), ("human", hum)):
    g = genome(path)
    for f in fracs:
        m = round(n * f / (100 - f))
        s1, s2 = sim(g, m, f"{donor}{f}")
        emit(f"{sid}__{donor}_{f}", h1 + s1, h2 + s2)
print("\n".join(rows))
