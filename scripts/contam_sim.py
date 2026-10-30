#!/usr/bin/env python3
"""Contamination spike-in inputs (#41): known fungal assemblies plus bacterial or human DNA, with a truth table.

usage: contam_sim.py HOST_DIR HOST_META.tsv DONOR_DIR OUT_DIR
  HOST_DIR/<accession>.fna for the accessions in HOST_META (stress_hosts.meta.tsv columns); DONOR_DIR holds bacteria.fna (E. coli +
  B. subtilis, concatenated) and human.fna (one human chromosome). Writes OUT_DIR/genomes/*.fa.gz and OUT_DIR/truth.tsv.
Variants per host: the host alone (frac 0) and host + donor with the donor a given share of the total bases, in 100 kb pieces.
Bacteria: 1, 2, 5, 10, 15% (two genomes are ~9 Mb, which is the most that is available without repeating sequence);
human: 1, 5, 10, 20, 30%. Truth = the host lineage; columns donor and frac say what was added.
"""
import csv, gzip, os, random, sys

hdir, meta, ddir, out = sys.argv[1:5]
os.makedirs(f"{out}/genomes", exist_ok=True)
rows = {r["id"]: r for r in csv.DictReader(open(meta), delimiter="\t")}
rng = random.Random(11)
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
FRACS = {"bacteria": (1, 2, 5, 10, 15), "human": (1, 5, 10, 20, 30)}


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


donors = {d: pieces(read(f"{ddir}/{d}.fna"), 100_000) for d in FRACS}
for d, p in donors.items():
    print(d, sum(map(len, p)) / 1e6, "Mb available", file=sys.stderr)
truth = csv.writer(open(f"{out}/truth.tsv", "w"), delimiter="\t", lineterminator="\n")
truth.writerow(["id", "host", "donor", "frac"] + RANKS)
n = 0
for h in [i for i in rows if os.path.exists(f"{hdir}/{i}.fna")]:
    hp = pieces(read(f"{hdir}/{h}.fna"), 100_000)
    total = sum(map(len, hp))
    lin = [rows[h][k] for k in RANKS]
    write(f"{h}__none_0", hp); truth.writerow([f"{h}__none_0", h, "none", 0] + lin); n += 1
    for d, dp in donors.items():
        dp = dp[:]; rng.shuffle(dp)
        for pct in FRACS[d]:
            need, acc, take = total * pct / (100 - pct), 0, []
            for p in dp:
                if acc >= need: break
                take.append(p); acc += len(p)
            if acc < need * 0.95:
                print(f"skip {h} {d} {pct}%: only {acc / 1e6:.1f} of {need / 1e6:.1f} Mb of donor", file=sys.stderr); continue
            name = f"{h}__{d}_{pct}"
            write(name, hp + take); truth.writerow([name, h, d, pct] + lin); n += 1
print(n, "inputs", file=sys.stderr)
