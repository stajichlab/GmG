#!/usr/bin/env python3
"""Pick a stratified non-fungal eukaryote decoy panel (#79) from `datasets summary genome taxon 2759 --reference`.

One genome per class (else phylum), phyla taken round-robin so the panel spreads across the tree, existing decoys first,
fungi (taxid 4751) excluded, genomes over --max-mb skipped to keep downloads small. Writes decoys.tsv rows.

  datasets summary genome taxon 2759 --reference --limit all --as-json-lines > euk.jsonl
  select_decoys.py euk.jsonl TAXDUMP_DIR data/decoys/decoys.tsv -n 90 > panel.tsv
"""
import argparse, json, sys, collections

ap = argparse.ArgumentParser()
ap.add_argument("jsonl"); ap.add_argument("taxdump"); ap.add_argument("current")
ap.add_argument("-n", type=int, default=90); ap.add_argument("--max-mb", type=float, default=1500)
ap.add_argument("--min-mb", type=float, default=5)
a = ap.parse_args()

parent, rank, name = {}, {}, {}
for l in open(f"{a.taxdump}/nodes.dmp"):
    f = [x.strip() for x in l.split("|")]
    parent[int(f[0])] = int(f[1]); rank[int(f[0])] = f[2]
for l in open(f"{a.taxdump}/names.dmp"):
    f = [x.strip() for x in l.split("|")]
    if f[3] == "scientific name": name[int(f[0])] = f[1]

def lineage(t):
    out = {}
    while t in parent and t != parent[t]:
        out.setdefault(rank[t], t); t = parent[t]
    return out

cur = [l.rstrip("\n").split("\t") for l in open(a.current)][1:]
have = {int(r[1]) for r in cur}
have_keys = set()
for r in cur:
    L = lineage(int(r[1])); have_keys.add(L.get("class") or L.get("phylum"))

cands = collections.defaultdict(list)
for l in open(a.jsonl):
    g = json.loads(l)
    acc = g["accession"]; t = int(g["organism"]["tax_id"])
    mb = int(g["assembly_stats"]["total_sequence_length"]) / 1e6
    L = lineage(t)
    if 4751 in L.values() or 2759 not in L.values():   # fungi out; Eukaryota is a "domain" in the current taxdump
        continue
    if not (a.min_mb <= mb <= a.max_mb) or not acc.startswith("GCF_"):
        continue
    key = L.get("class") or L.get("phylum") or L.get("kingdom") or L.get("clade")
    if key is None or key in have_keys: continue
    ph = L.get("phylum") or L.get("kingdom") or key
    cands[(ph, key)].append((mb, acc, t, g["organism"]["organism_name"], L))

byph = collections.defaultdict(list)
for (ph, key), v in cands.items():
    v.sort(); byph[ph].append((key, v[0]))   # smallest eligible genome: cheapest download
picked = []
while len(picked) < a.n - len(cur) and any(byph.values()):
    for ph in sorted(byph, key=lambda p: -len(byph[p])):
        if byph[ph] and len(picked) < a.n - len(cur):
            picked.append(byph[ph].pop(0))
print("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies")
for r in cur: print("\t".join(r))
tot = 0
for key, (mb, acc, t, nm, L) in picked:
    tot += mb
    print(f"{acc}\t{t}\t\t\t\t\t\t{nm}")
print(f"# {len(cur)} existing + {len(picked)} new, {tot/1000:.1f} Gb of new genomes, {len(byph)} phyla", file=sys.stderr)
