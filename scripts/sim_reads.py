#!/usr/bin/env python3
"""Simulate single-end Illumina-like reads from a genome FASTA (substitution errors only).

usage: sim_reads.py GENOME.fna OUT.fastq.gz COVERAGE [--len 150] [--err 0.005] [--seed 1]
Reads start uniformly over contigs (prob. proportional to length), random strand.
"""
import sys, gzip, argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("genome"); ap.add_argument("out"); ap.add_argument("coverage", type=float)
ap.add_argument("--len", type=int, default=150); ap.add_argument("--err", type=float, default=0.005)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)

seqs, cur = [], []
for line in open(a.genome):
    if line.startswith(">"):
        if cur: seqs.append("".join(cur)); cur = []
    else:
        cur.append(line.strip().upper())
if cur: seqs.append("".join(cur))
seqs = [s for s in seqs if len(s) >= a.len]
lens = np.array([len(s) for s in seqs], dtype=np.int64)
total = int(lens.sum())
n_reads = int(a.coverage * total / a.len)
tbl = np.full(256, 78, dtype=np.uint8)  # N
for i, b in enumerate(b"ACGT"): tbl[b] = i
comp = np.array([3, 2, 1, 0], dtype=np.uint8)
letters = np.frombuffer(b"ACGTN", dtype=np.uint8)

contig = rng.choice(len(seqs), size=n_reads, p=lens / total)
out = gzip.open(a.out, "wt", compresslevel=3)
qual = "I" * a.len
rid = 0
for ci in np.unique(contig):
    n = int((contig == ci).sum())
    arr = tbl[np.frombuffer(seqs[ci].encode(), dtype=np.uint8)]
    starts = rng.integers(0, lens[ci] - a.len + 1, size=n)
    idx = starts[:, None] + np.arange(a.len)[None, :]
    reads = arr[idx]                                   # n x len, codes 0-3 (78 for N)
    rc = rng.random(n) < 0.5
    flip = reads[rc][:, ::-1]
    ok = flip < 4
    flip = np.where(ok, comp[np.minimum(flip, 3)], flip)
    reads[rc] = flip
    err = rng.random(reads.shape) < a.err
    sub = (reads + rng.integers(1, 4, size=reads.shape).astype(np.uint8)) % 4
    reads = np.where(err & (reads < 4), sub, reads)
    txt = np.where(reads < 4, letters[np.minimum(reads, 3)], letters[4])
    for r in txt:
        out.write(f"@sim{rid}\n{r.tobytes().decode()}\n+\n{qual}\n"); rid += 1
out.close()
print(f"{a.genome}: {total/1e6:.1f} Mb, {n_reads} reads, {a.coverage}x", file=sys.stderr)
