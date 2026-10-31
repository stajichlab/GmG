#!/usr/bin/env python3
"""Simulate ONT/PacBio-like long reads (substitutions + insertions + deletions) from a genome FASTA, at several depths in one pass (numpy only).

usage: sim_long_reads.py GENOME.fa[.gz] OUT_PREFIX --err 0.05 [--depths 0.25,0.5,1,2,5,10,30] [--sub 0.4 --ins 0.3 --del 0.3]
                         [--median-len 12000 --sigma 0.5 --min-len 2000 --max-len 60000] [--seed 1]
Writes OUT_PREFIX_d<depth>.fastq for every depth. The files are nested prefixes of one random read stream (the 0.25x file is the first reads of the
30x file), so a lower depth is a random subsample of a higher one. Reads: log-normal length (median 12 kb, sigma 0.5: the central 95% are 4.5-32 kb, mean 13.6 kb), uniform start over contigs (probability proportional to length), random strand. Errors: every position of the read is an
error with probability --err, of type substitution / insertion (one random base) / deletion in the proportions --sub/--ins/--del (default 0.4/0.3/0.3: ONT-like;
HiFi is mostly substitutions). Uniform and independent along the read: real errors are context dependent (homopolymers) and clustered, so the k-mer survival
here, (1 - err)^21 for substitutions, is a lower bound of the real one. Qualities are a constant 'I'. This is a small stand-in for PBSIM3/Badread, which are not
installed on the HPCC; scripts/long_reads_measure.py (minimap2 alignment) checks the realised error rate.
"""
import sys, argparse, gzip
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("genome"); ap.add_argument("prefix")
ap.add_argument("--err", type=float, required=True)
ap.add_argument("--depths", default="0.25,0.5,1,2,5,10,30")
ap.add_argument("--sub", type=float, default=0.4); ap.add_argument("--ins", type=float, default=0.3); ap.add_argument("--del", dest="dele", type=float, default=0.3)
ap.add_argument("--median-len", type=float, default=12000); ap.add_argument("--sigma", type=float, default=0.5)
ap.add_argument("--min-len", type=int, default=2000); ap.add_argument("--max-len", type=int, default=60000)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)

tbl = np.full(256, 255, dtype=np.uint8)
for i, b in enumerate(b"ACGT"): tbl[b] = i; tbl[b + 32] = i
op = gzip.open if a.genome.endswith(".gz") else open
seqs, cur = [], []
with op(a.genome, "rt") as f:
    for line in f:
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip())
if cur: seqs.append("".join(cur))
contigs = []
for s in seqs:
    if len(s) < a.min_len: continue
    c = tbl[np.frombuffer(s.encode(), dtype=np.uint8)]
    bad = c == 255
    c[bad] = rng.integers(0, 4, size=int(bad.sum()))          # N and IUPAC codes: a random base
    contigs.append(c)
lens = np.array([len(c) for c in contigs], dtype=np.int64)
total = int(lens.sum())
depths = [float(x) for x in a.depths.split(",")]
targets = [d * total for d in depths]
outs = [open(f"{a.prefix}_d{d:g}.fastq", "w") for d in depths]
done = [False] * len(depths)
cum = 0
letters = np.frombuffer(b"ACGT", dtype=np.uint8)
comp = np.array([3, 2, 1, 0], dtype=np.uint8)
ps, pi, pd = np.array([a.sub, a.ins, a.dele]) / (a.sub + a.ins + a.dele)
n = 0
tot_err = tot_len = 0
while not all(done):
    L = int(np.clip(rng.lognormal(np.log(a.median_len), a.sigma), a.min_len, a.max_len))
    ci = int(rng.choice(len(contigs), p=lens / total))
    c = contigs[ci]
    if L > len(c): L = len(c)
    st = int(rng.integers(0, len(c) - L + 1))
    t = c[st:st + L].copy()
    if rng.random() < 0.5:
        t = comp[t[::-1]]
    ev = rng.random(L) < a.err
    u = rng.random(L)
    sub = ev & (u < ps); ins = ev & (u >= ps) & (u < ps + pi); dele = ev & (u >= ps + pi)
    t[sub] = (t[sub] + rng.integers(1, 4, size=int(sub.sum())).astype(np.uint8)) % 4
    counts = np.ones(L, dtype=np.int64); counts[dele] = 0; counts[ins] = 2
    out = np.repeat(t, counts)
    if ins.any():
        pos = np.cumsum(counts)[ins] - 1                       # the second copy of an inserted position gets a random base
        out[pos] = rng.integers(0, 4, size=len(pos)).astype(np.uint8)
    seq = letters[out].tobytes().decode()
    rec = f"@sim{n}_c{ci}_{st}\n{seq}\n+\n{'I' * len(seq)}\n"
    n += 1; cum += len(seq); tot_err += int(ev.sum()); tot_len += L
    for k, (tg, fh) in enumerate(zip(targets, outs)):
        if not done[k]:
            fh.write(rec)
            if cum >= tg: done[k] = True; fh.close()
print(f"{a.genome}: {total/1e6:.2f} Mb, err {a.err}, {n} reads, mean read {cum/n:.0f}, realised event rate {tot_err/tot_len:.4f}, depths {a.depths}", file=sys.stderr)
