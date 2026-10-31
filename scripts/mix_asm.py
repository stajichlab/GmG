#!/usr/bin/env python3
"""Two-fungus assembly mixtures for the mixture detector (#31).

usage: mix_asm.py TRUTH_ALL.tsv GENOME_DIR OUT_DIR [--per-class N] [--seed S]
  TRUTH_ALL.tsv  the RefSeq truth set (GmG_runs/truthset/truth.all.tsv: lineage, bfd_id, in_bfd, total_len)
  GENOME_DIR     <id>.fna.gz of every truth-set genome
Picks N genome pairs (default 12) per class from truth-set genomes that are in the reference and 8-60 Mb, every genome used once:
  same_genus (two species of a genus), same_order (two genera of an order), diff_phylum (two phyla).
For each pair A (majority) and B (minority) writes assemblies in which B makes up 0.1, 0.5, 1, 2, 5, 10, 20, 30 and 50% of the bases (both
genomes in 100 kb pieces, as scripts/stress2_sim.py), plus A alone and B alone as pure controls.
Writes OUT_DIR/genomes/*.fa.gz, OUT_DIR/sheet.tsv (id, path, expected), OUT_DIR/truth.tsv and OUT_DIR/exclude_ids.txt: the reference
ids of every parent, to be hidden with `gmg run --exclude-ids` (the parents are in the reference; their relatives stay).
"""
import argparse, csv, gzip, os, random

ap = argparse.ArgumentParser()
ap.add_argument("truth"); ap.add_argument("gdir"); ap.add_argument("out")
ap.add_argument("--per-class", type=int, default=12)
ap.add_argument("--seed", type=int, default=31)
a = ap.parse_args()
rng = random.Random(a.seed)
L = ["phylum", "class", "order", "family", "genus", "species"]
SHARES = [0.1, 0.5, 1, 2, 5, 10, 20, 30, 50]
os.makedirs(f"{a.out}/genomes", exist_ok=True)

rows = [r for r in csv.DictReader(open(a.truth), delimiter="\t")
        if r["in_bfd"] == "yes" and r["bfd_id"] and r["species"] and 8e6 <= float(r["total_len"]) <= 60e6
        and os.path.exists(f"{a.gdir}/{r['id']}.fna.gz")]
rng.shuffle(rows)
used = set()


def find(kind):
    for x in rows:
        if x["id"] in used: continue
        for y in rows:
            if y["id"] in used or y["id"] == x["id"]: continue
            if kind == "same_genus" and x["genus"] and x["genus"] == y["genus"] and x["species"] != y["species"]: return x, y
            if kind == "same_order" and x["order"] and x["order"] == y["order"] and x["genus"] != y["genus"] and x["genus"] and y["genus"]: return x, y
            if kind == "diff_phylum" and x["phylum"] != y["phylum"]: return x, y
    return None


def read(path):
    seqs, cur = [], []
    for line in gzip.open(path, "rt"):
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
        if acc + len(p) > want: p = p[:int(want - acc)]
        if not p: break
        keep.append(p); acc += len(p)
    return keep


def write(name, contigs):
    with gzip.open(f"{a.out}/genomes/{name}.fa.gz", "wt", compresslevel=1) as f:
        for i, c in enumerate(contigs):
            f.write(f">{name}_{i}\n")
            for j in range(0, len(c), 80): f.write(c[j:j + 80] + "\n")


sheet = csv.writer(open(f"{a.out}/sheet.tsv", "w"), delimiter="\t", lineterminator="\n")
sheet.writerow(["id", "path", "expected"])
truth = csv.writer(open(f"{a.out}/truth.tsv", "w"), delimiter="\t", lineterminator="\n")
truth.writerow(["id", "kind", "pair", "b_pct", "a", "b"] + L + ["b_" + x for x in L])
excl = []
for kind in ("same_genus", "same_order", "diff_phylum"):
    for k in range(a.per_class):
        p = find(kind)
        if not p: print("no more pairs for", kind); break
        A, B = p
        if rng.random() < 0.5: A, B = B, A
        used.update([A["id"], B["id"]]); excl += [A["bfd_id"], B["bfd_id"]]
        pa, pb = pieces(read(f"{a.gdir}/{A['id']}.fna.gz")), pieces(read(f"{a.gdir}/{B['id']}.fna.gz"))
        la, lb = sum(map(len, pa)), sum(map(len, pb))
        lin = [A[x] for x in L] + [B[x] for x in L]
        pid = f"{kind}{k:02d}"
        for pct in [0] + SHARES + [100]:
            q = pct / 100
            if pct == 0: ctg = pa
            elif pct == 100: ctg = pb
            else:
                total = min(la / (1 - q), lb / q)
                ctg = pick(pa, total * (1 - q)) + pick(pb, total * q)
            name = f"{pid}__B{pct:g}"
            write(name, ctg)
            sheet.writerow([name, f"{a.out}/genomes/{name}.fa.gz", ""])
            truth.writerow([name, kind, pid, pct, A["id"], B["id"]] + lin)
open(f"{a.out}/exclude_ids.txt", "w").write("\n".join(excl) + "\n")
print(len(used) // 2, "pairs")
