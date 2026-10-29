#!/usr/bin/env python3
"""Stress-test inputs (issue #41): degrade known genomes and write a truth table.

usage: stress_sim.py GENOME_DIR META.tsv OUT_DIR
  GENOME_DIR holds <accession>.fna for accessions in META (post_freeze.meta.tsv columns). Writes OUT_DIR/genomes/*.fa.gz and OUT_DIR/truth.tsv
Variants per host genome (truth = the host lineage):
  complete{10,25,50,75}  random 100 kb pieces totalling that share of the genome (completeness gradient)
  frag5k                 whole genome cut into 5 kb contigs (N50 flag test; sketch content unchanged)
  contam{5,20,50}        host plus donor bases from a genome of the other phylum, donor = that share of the total
Negative controls: random_acgt (uniform random DNA), repeat (low-complexity tandem repeat); truth phylum blank.
"""
import csv, gzip, os, random, sys

gdir, meta, out = sys.argv[1:4]
os.makedirs(f"{out}/genomes", exist_ok=True)
rows = {r["id"]: r for r in csv.DictReader(open(meta), delimiter="\t")}
rng = random.Random(7)


def read(path):
    seqs, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else:
            cur.append(line.strip().upper())
    if cur: seqs.append("".join(cur))
    return seqs


def pieces(seqs, size):
    return [s[i:i + size] for s in seqs for i in range(0, len(s), size) if len(s[i:i + size]) >= size // 2]


def write(name, contigs):
    with gzip.open(f"{out}/genomes/{name}.fa.gz", "wt", compresslevel=1) as f:
        for i, c in enumerate(contigs):
            f.write(f">{name}_{i}\n")
            for j in range(0, len(c), 80): f.write(c[j:j + 80] + "\n")


ids = [i for i in rows if os.path.exists(f"{gdir}/{i}.fna")]
seq = {i: read(f"{gdir}/{i}.fna") for i in ids}
truth = csv.writer(open(f"{out}/truth.tsv", "w"), delimiter="\t")
truth.writerow(["id", "host", "variant", "param", "phylum", "class", "order", "family", "genus", "species"])


def emit(name, host, var, param, contigs):
    write(name, contigs)
    r = rows[host] if host else {}
    truth.writerow([name, host, var, param] + [r.get(k, "") for k in ["phylum", "class", "order", "family", "genus", "species"]])


for h in ids:
    pcs = pieces(seq[h], 100_000)
    total = sum(len(p) for p in pcs)
    for pct in (10, 25, 50, 75):
        rng.shuffle(pcs)
        keep, acc = [], 0
        for p in pcs:
            if acc >= total * pct / 100: break
            keep.append(p); acc += len(p)
        emit(f"{h}__complete{pct}", h, "complete", pct, keep)
    emit(f"{h}__frag5k", h, "frag5k", 5000, pieces(seq[h], 5000))
    donors = [d for d in ids if rows[d]["phylum"] != rows[h]["phylum"]]
    d = rng.choice(donors)
    dp = pieces(seq[d], 100_000); rng.shuffle(dp)
    for pct in (5, 20, 50):
        need, acc, take = total * pct / (100 - pct), 0, []
        for p in dp:
            if acc >= need: break
            take.append(p); acc += len(p)
        emit(f"{h}__contam{pct}", h, "contam", pct, pieces(seq[h], 100_000) + take)
emit("neg_random_acgt", "", "negative", "random", ["".join(rng.choice("ACGT") for _ in range(100_000)) for _ in range(150)])
emit("neg_repeat", "", "negative", "repeat", ["ACGTTGCAAT" * 10_000] * 150)
print(len(ids), "hosts;", sum(1 for _ in open(f"{out}/truth.tsv")) - 1, "inputs", file=sys.stderr)
