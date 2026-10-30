#!/usr/bin/env python3
"""Chimeras and hybrids (#41, #45): usage stress2_sim.py PAIRS.tsv GENOME_DIR OUT_DIR

PAIRS.tsv from stress2_pairs.py; GENOME_DIR/<accession>.fna. Writes OUT_DIR/genomes/*.fa.gz and OUT_DIR/truth.tsv.
Per pair (A, B):
  chimera 50:50, 70:30, 90:10  A is the majority; both genomes in 100 kb pieces, the amounts chosen so the shares are exact
                               and as much of the smaller genome as possible is used (total = min(|A|/p, |B|/(1-p)))
  hybrid (same_genus pairs)    A and B in full: what an allodiploid / interspecies hybrid assembly contains
Truth = the lineage of A (the majority); b_* columns hold the minority; at 50:50 either parent is a correct call.
"""
import csv, gzip, os, random, sys

pairs_f, gdir, out = sys.argv[1:4]
os.makedirs(f"{out}/genomes", exist_ok=True)
rng = random.Random(45)
L = ["phylum", "class", "order", "family", "genus", "species"]


def read(path):
    seqs, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip().upper())
    if cur: seqs.append("".join(cur))
    return seqs


def pieces(seqs, size=100_000):
    return [s[i:i + size] for s in seqs for i in range(0, len(s), size) if len(s[i:i + size]) >= size // 2]


def pick(pcs, want):
    pcs = pcs[:]; rng.shuffle(pcs)
    keep, acc = [], 0
    for p in pcs:
        if acc >= want: break
        if acc + len(p) > want:
            p = p[:int(want - acc)]   # the last piece is cut so the share is exact
        if not p:                     # a remainder below one base: an empty record would not parse
            break
        keep.append(p); acc += len(p)
    return keep


def write(name, contigs):
    with gzip.open(f"{out}/genomes/{name}.fa.gz", "wt", compresslevel=1) as f:
        for i, c in enumerate(contigs):
            f.write(f">{name}_{i}\n")
            for j in range(0, len(c), 80): f.write(c[j:j + 80] + "\n")


truth = csv.writer(open(f"{out}/truth.tsv", "w"), delimiter="\t", lineterminator="\n")
truth.writerow(["id", "kind", "variant", "pct", "a", "b"] + L + ["b_" + x for x in L])
n = 0
for r in csv.DictReader(open(pairs_f), delimiter="\t"):
    a, b = r["a"], r["b"]
    if not (os.path.exists(f"{gdir}/{a}.fna") and os.path.exists(f"{gdir}/{b}.fna")):
        print("missing genome for", a, b, file=sys.stderr); continue
    pa, pb = pieces(read(f"{gdir}/{a}.fna")), pieces(read(f"{gdir}/{b}.fna"))
    la, lb = sum(map(len, pa)), sum(map(len, pb))
    lin = [r["a_" + x] for x in L] + [r["b_" + x] for x in L]
    for pct in (50, 70, 90):
        p = pct / 100
        total = min(la / p, lb / (1 - p))
        name = f"{a}__{b}__{r['kind']}_chim{pct}"
        write(name, pick(pa, total * p) + pick(pb, total * (1 - p)))
        truth.writerow([name, r["kind"], "chimera", pct, a, b] + lin); n += 1
    if r["kind"] == "same_genus":
        name = f"{a}__{b}__{r['kind']}_hybrid"
        write(name, pa + pb)
        truth.writerow([name, r["kind"], "hybrid", 100, a, b] + lin); n += 1
print(n, "inputs", file=sys.stderr)
