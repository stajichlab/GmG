#!/usr/bin/env python3
"""All NCBI fungal assemblies released after the BFD freeze that are not in BFD, at every quality level (no N50/size filter).

usage: select_post_freeze.py DATA_DIR OUT_META.tsv [RELEASED_AFTER=2025-05-01]
  DATA_DIR has fungi_all.tsv (datasets + dataformat), bfd_meta.tsv and rankedlineage.dmp (see scripts/select_external.py).
Truth lineage comes from the NCBI taxonomy (rankedlineage.dmp). Columns species_in_bfd / genus_in_bfd tell known from novel taxa;
level, n50, total_len, release allow accuracy to be broken down by assembly quality.
"""
import sys, collections

d, out = sys.argv[1], sys.argv[2]
after = sys.argv[3] if len(sys.argv) > 3 else "2025-05-01"
bfd = [l.split("\t") for l in open(f"{d}/bfd_meta.tsv").read().splitlines()[1:]]
bfd_nums = {r[0].split("_")[1].split(".")[0] for r in bfd}
bfd_species = {r[7] for r in bfd if r[7]}
bfd_genera = {r[6] for r in bfd if r[6]}
lin = {}
with open(f"{d}/rankedlineage.dmp") as fh:
    for line in fh:
        f = [x.strip() for x in line.rstrip("\n").split("\t|")]
        lin[f[0]] = dict(species=f[2], genus=f[3], family=f[4], order=f[5], cls=f[6], phylum=f[7])
rows = [l.split("\t") for l in open(f"{d}/fungi_all.tsv").read().splitlines()[1:]]
keep, skipped = [], collections.Counter()
for r in rows:
    acc, name, taxid, rel, level, tlen, n50 = r[:7]
    if acc.split("_")[1].split(".")[0] in bfd_nums:
        skipped["in BFD"] += 1; continue
    if rel < after:
        skipped["released before cutoff"] += 1; continue
    l = lin.get(taxid)
    if not l or not l["phylum"]:
        skipped["no phylum in NCBI lineage"] += 1; continue
    keep.append((acc, taxid, l, rel, level, n50 or "0", tlen or "0"))
keep.sort()
with open(out, "w") as m:
    m.write("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\tspecies_in_bfd\tgenus_in_bfd\trelease\tlevel\tn50\ttotal_len\n")
    for acc, taxid, l, rel, level, n50, tlen in keep:
        m.write("\t".join([acc, taxid, l["phylum"], l["cls"], l["order"], l["family"], l["genus"], l["species"],
                           str(l["species"] in bfd_species and l["species"] != ""), str(l["genus"] in bfd_genera and l["genus"] != ""),
                           rel, level, n50, tlen]) + "\n")
print(f"selected {len(keep)} assemblies; skipped {dict(skipped)}")
c = collections.Counter(k[2]["phylum"] for k in keep)
print("phyla:", dict(c.most_common(8)))
print("species named:", sum(1 for k in keep if k[2]["species"]), "| species in BFD:", sum(1 for k in keep if k[2]["species"] in bfd_species and k[2]["species"]),
      "| genus not in BFD:", sum(1 for k in keep if k[2]["genus"] and k[2]["genus"] not in bfd_genera))
print("levels:", dict(collections.Counter(k[4] for k in keep)))
