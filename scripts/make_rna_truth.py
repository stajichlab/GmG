#!/usr/bin/env python3
"""Stratified sample of BFD RNA-seq read sets with species truth, for the RNA validation (issue #43).

usage: make_rna_truth.py SAMPLES_CSV RNASEQ_READS_DIR OUT_DIR [N=300] [MAX_PER_GENUS=3]
RNASEQ_READS_DIR holds <Species_name>_norm_{R1,R2,SE}.fastq.gz (reads pooled from up to 5 SRA runs, cleaned and normalised for assembly).
Truth lineage = the BFD samples.csv lineage of that species. Takes R1 if present, else SE. Stratified by phylum (proportional, at least
MIN_PER_PHYLUM=8 where available), at most MAX_PER_GENUS per genus, fixed seed.
writes OUT_DIR/rna.truth.tsv and OUT_DIR/rna.list (paths, one per line).
"""
import csv, os, random, sys, collections

samples, rdir, out = sys.argv[1:4]
n = int(sys.argv[4]) if len(sys.argv) > 4 else 300
per_genus = int(sys.argv[5]) if len(sys.argv) > 5 else 3
files = set(os.listdir(rdir))
lin = {}
with open(samples, newline="") as fh:
    for r in csv.DictReader(fh):
        sp = r["SPECIES"].strip()
        if sp: lin.setdefault(sp.replace(" ", "_"), r)
cands = []
for f in sorted(files):
    if not f.endswith("_norm_R1.fastq.gz"): continue
    sp = f[: -len("_norm_R1.fastq.gz")]
    if sp in lin: cands.append((sp, f, lin[sp]))
for f in sorted(files):   # single-end only species
    if f.endswith("_norm_SE.fastq.gz"):
        sp = f[: -len("_norm_SE.fastq.gz")]
        if sp in lin and not any(c[0] == sp for c in cands): cands.append((sp, f, lin[sp]))
random.seed(7); random.shuffle(cands)
by_ph = collections.defaultdict(list)
for c in cands: by_ph[c[2]["PHYLUM"] or "(none)"].append(c)
tot = sum(len(v) for v in by_ph.values())
pick, genus_n = [], collections.Counter()
for ph, lst in by_ph.items():
    q = max(min(8, len(lst)), round(n * len(lst) / tot))
    taken = 0
    for sp, f, r in lst:
        if taken >= q: break
        if genus_n[r["GENUS"]] >= per_genus: continue
        genus_n[r["GENUS"]] += 1; taken += 1; pick.append((sp, f, r))
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "rna.truth.tsv"), "w") as t, open(os.path.join(out, "rna.list"), "w") as l:
    t.write("id\tphylum\tclass\torder\tfamily\tgenus\tspecies\tfile\n")
    for sp, f, r in sorted(pick):
        rid = f.replace(".fastq.gz", "")
        t.write("\t".join([rid, r["PHYLUM"], r["CLASS"], r["ORDER"], r["FAMILY"], r["GENUS"], r["SPECIES"], f]) + "\n")
        l.write(os.path.join(rdir, f) + "\n")
print(len(pick), "read sets;", dict(collections.Counter(r["PHYLUM"] for _, _, r in pick).most_common(6)))
