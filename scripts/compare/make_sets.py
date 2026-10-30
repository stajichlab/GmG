#!/usr/bin/env python3
"""Query sets for the tool comparison (#46). usage: make_sets.py REPO_DIR OUT_DIR SEED
Writes OUT_DIR/s1_post_freeze.tsv (post-freeze NCBI assemblies, independent of the reference: all with the species in the reference,
all with the genus not in it, a random sample of the rest, ~300) and OUT_DIR/s2_truth.tsv (RefSeq truth-set genomes that are in the index,
~300: all outside Ascomycota/Basidiomycota, 128 of each of those two; the query is the BFD copy of the genome, scored leave-one/species/
genus-out by dropping hits after the search). Stdlib only."""
import csv, os, random, sys

repo, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
rng = random.Random(seed)
os.makedirs(out, exist_ok=True)
R = ["phylum", "class", "order", "family", "genus", "species"]
pf = list(csv.DictReader(open(f"{repo}/data/validation/post_freeze.meta.tsv"), delimiter="\t"))
sp_in = [r for r in pf if r["species_in_bfd"] == "True"]
gn_out = [r for r in pf if r["genus_in_bfd"] == "False"]
rest = [r for r in pf if r not in sp_in and r not in gn_out]
pick = {r["id"]: r for r in sp_in + gn_out}
for r in rng.sample(rest, min(len(rest), 300 - len(pick))):
    pick[r["id"]] = r
cols = ["id"] + R + ["species_in_bfd", "genus_in_bfd", "level", "n50", "total_len"]
with open(f"{out}/s1_post_freeze.tsv", "w") as o:
    w = csv.writer(o, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for r in sorted(pick.values(), key=lambda r: r["id"]):
        w.writerow([r[c] for c in cols])
tr = [r for r in csv.DictReader(open(f"{repo}/data/validation/truth_set.list.tsv"), delimiter="\t") if r["in_bfd"] == "yes" and r["bfd_id"]]
big = [r for r in tr if r["phylum"] in ("Ascomycota", "Basidiomycota")]
small = [r for r in tr if r not in big]
sel = small + rng.sample([r for r in big if r["phylum"] == "Ascomycota"], 128) + rng.sample([r for r in big if r["phylum"] == "Basidiomycota"], 128)
cols2 = ["id", "bfd_id"] + R
with open(f"{out}/s2_truth.tsv", "w") as o:
    w = csv.writer(o, delimiter="\t", lineterminator="\n")
    w.writerow(cols2)
    for r in sorted(sel, key=lambda r: r["id"]):
        w.writerow([r[c] for c in cols2])
print(len(pick), "post-freeze;", len(sel), "truth-set genomes", file=sys.stderr)
