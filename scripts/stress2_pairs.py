#!/usr/bin/env python3
"""Pick genome pairs for the chimera and hybrid stress tests (#41, #45): usage stress2_pairs.py POST_FREEZE.meta.tsv OUT.tsv [N_PER_KIND]

Chromosome/complete assemblies of 8-60 Mb with a named genus, genus in the reference. Kinds (the relation of the two genomes):
same_genus (different species), same_order (different genus), same_class (different order). Seeded, one pair per genus pair where possible.
"""
import csv, random, sys, collections

meta, out = sys.argv[1:3]
n_each = int(sys.argv[3]) if len(sys.argv) > 3 else 8
rng = random.Random(41)
rows = [r for r in csv.DictReader(open(meta), delimiter="\t")
        if r["level"] in ("Chromosome", "Complete Genome", "Complete") and r["genus"] and r["genus_in_bfd"] == "True"
        and 8e6 <= float(r["total_len"] or 0) <= 60e6]
rng.shuffle(rows)
by = collections.defaultdict(list)
for r in rows:
    by[r["class"]].append(r)


def species_key(r):
    return r["species"] or r["taxid"]


pairs, used = [], set()


def take(kind, ok, same):
    got = 0
    cands = [(a, b) for a in rows for b in rows if a["id"] < b["id"] and a["phylum"] == b["phylum"] and ok(a, b)]
    rng.shuffle(cands)
    seen = set()
    for a, b in cands:
        key = tuple(sorted((a["genus"], b["genus"]))) if kind != "same_genus" else (a["genus"],)
        if key in seen or a["id"] in used or b["id"] in used:
            continue
        seen.add(key); used.update((a["id"], b["id"]))
        pairs.append((kind, a, b)); got += 1
        if got >= n_each: break
    print(kind, got, file=sys.stderr)


take("same_genus", lambda a, b: a["genus"] == b["genus"] and species_key(a) != species_key(b), True)
take("same_order", lambda a, b: a["order"] == b["order"] and a["genus"] != b["genus"], False)
take("same_class", lambda a, b: a["class"] == b["class"] and a["order"] != b["order"], False)
L = ["phylum", "class", "order", "family", "genus", "species"]
with open(out, "w") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["kind", "a", "b", "a_len", "b_len"] + ["a_" + x for x in L] + ["b_" + x for x in L])
    for kind, a, b in pairs:
        w.writerow([kind, a["id"], b["id"], a["total_len"], b["total_len"]] + [a[x] for x in L] + [b[x] for x in L])
