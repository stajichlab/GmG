#!/usr/bin/env python3
"""Truth table from lab labels via the NCBI taxdump (used by the reads validations, #42/#45).

usage: name_truth.py LABELS.tsv TAXDUMP_DIR OUT.tsv [--index-meta META.tsv]     LABELS: columns sample, name ("Genus species ...")
Lineage (phylum..genus) = that of the name's genus, found in names.dmp as a scientific name, synonym, includes or equivalent name, so
old genera (Absidia, Hypocrea) resolve to the current one. Species is left blank: lab species labels are not curated. Writes id, phylum, class,
order, family, genus, species (blank), label_genus, and with --index-meta genus_in_ref (True when the resolved genus is in the index table).
Labels whose genus is not in the taxdump are dropped (counted on stderr).
"""
import sys, re, csv, collections

args = [a for a in sys.argv[1:] if not a.startswith("--")]
meta = sys.argv[sys.argv.index("--index-meta") + 1] if "--index-meta" in sys.argv else None
if meta: args.remove(meta)
lab_p, td, out = args[:3]
RANKS = ["phylum", "class", "order", "family", "genus"]
labels = {r["sample"]: r["name"] for r in csv.DictReader(open(lab_p), delimiter="\t")}
genus_of = {}
for s, n in labels.items():
    n = re.sub(r"[\[\]]|\bcf\.?\b|\baff\.?\b", "", n).strip()
    w = n.split()[0].split("/")[0] if n else ""
    if w[:1].isupper(): genus_of[s] = w
want = set(genus_of.values())
by_name = collections.defaultdict(set)
for l in open(f"{td}/names.dmp"):
    f = [x.strip() for x in l.split("|")]
    if f[1] in want and f[3] in ("scientific name", "synonym", "includes", "equivalent name"): by_name[f[1]].add(int(f[0]))
parent, rank = {}, {}
for l in open(f"{td}/nodes.dmp"):
    f = [x.strip() for x in l.split("|")]
    parent[int(f[0])] = int(f[1]); rank[int(f[0])] = f[2]
anc = {}
for t in {t for ts in by_name.values() for t in ts}:
    x, d = t, {}
    while x in parent and parent[x] != x:
        d.setdefault(rank[x], x); x = parent[x]
    anc[t] = d
sci = {}
allt = {v for d in anc.values() for v in d.values()}
for l in open(f"{td}/names.dmp"):
    f = [x.strip() for x in l.split("|")]
    if f[3] == "scientific name" and int(f[0]) in allt: sci[int(f[0])] = f[1]
ref_genera = set()
if meta:
    ref_genera = {r["genus"] for r in csv.DictReader(open(meta), delimiter="\t") if r.get("genus")}


def truth(g):
    best = None
    for t in by_name.get(g, ()):
        d = anc[t]
        if rank[t] == "genus" or "genus" in d:
            lin = {r: sci.get(d.get(r), "") for r in RANKS}
            if best is None or len(lin["phylum"]) > len(best["phylum"]): best = lin
    return best


cols = ["id"] + RANKS + ["species", "label_genus"] + (["genus_in_ref"] if meta else [])
dropped = 0
with open(out, "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(cols)
    for s in labels:
        g = genus_of.get(s); tr = truth(g) if g else None
        if not tr or not tr["phylum"]: dropped += 1; continue
        w.writerow([s] + [tr[r] for r in RANKS] + ["", g] + ([str(tr["genus"] in ref_genera)] if meta else []))
print(len(labels) - dropped, "labels resolved;", dropped, "dropped (genus not in the taxdump)", file=sys.stderr)
