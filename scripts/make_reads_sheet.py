#!/usr/bin/env python3
"""Sample sheet for the DNA-reads titration (#42) from a folder of paired Illumina FASTQ named <Genus_species_STRAIN>_R1/_R2.f(ast)q.gz.

usage: make_reads_sheet.py INPUT_DIR OUT_DIR [PER_GENUS=2] [SEED=5]
Takes up to PER_GENUS pairs per genus (both files present, readable, >= 30 MB), skips broken links and non-FASTQ files, and writes
OUT_DIR/sheet.tsv (id, path, expected: two rows per sample, R1 and R2, merged into one sketch by `gmg run --sample-sheet`) and
OUT_DIR/labels.tsv (sample, name = the first two words of the file name: the lab's label, to be resolved by name_truth.py).
Read-only on INPUT_DIR.
"""
import os, re, random, sys, collections

src, out = sys.argv[1:3]
per = int(sys.argv[3]) if len(sys.argv) > 3 else 2
random.seed(int(sys.argv[4]) if len(sys.argv) > 4 else 5)
pat = re.compile(r"^(.*)_R([12])\.f(?:ast)?q\.gz$")
pairs = collections.defaultdict(dict)
for f in sorted(os.listdir(src)):
    m = pat.match(f)
    if m: pairs[m.group(1)][m.group(2)] = f
by_genus = collections.defaultdict(list)
for stem, d in pairs.items():
    if len(d) < 2: continue
    try:
        if min(os.path.getsize(os.path.join(src, d[k])) for k in "12") < 30_000_000: continue
    except OSError:
        continue   # broken link
    g = stem.split("_")[0]
    if g[:1].isupper(): by_genus[g].append(stem)
os.makedirs(out, exist_ok=True)
n = 0
with open(f"{out}/sheet.tsv", "w") as sh, open(f"{out}/labels.tsv", "w") as lb:
    sh.write("id\tpath\texpected\n"); lb.write("sample\tname\n")
    for g in sorted(by_genus):
        random.shuffle(by_genus[g])
        for stem in sorted(by_genus[g][:per]):
            sid = re.sub(r"[^A-Za-z0-9_.-]", "_", stem)
            name = " ".join(stem.split("_")[:2])
            for k in "12":
                sh.write(f"{sid}\t{os.path.join(src, pairs[stem][k])}\t\n")
            lb.write(f"{sid}\t{name}\n"); n += 1
print(n, "samples from", len(by_genus), "genera", file=sys.stderr)
