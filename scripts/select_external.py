#!/usr/bin/env python3
"""Pick fungal assemblies that are NOT in BFD (by accession number) for held-out testing.

usage: select_external.py DATA_DIR N_PER_PHYLUM_JSON OUT_PREFIX
  DATA_DIR has fungi_all.tsv (datasets+dataformat), bfd_meta.tsv, rankedlineage.dmp
  Lineage comes from the NCBI taxonomy (rankedlineage.dmp), like BFD's own columns.
Quality filter: N50 >= 20 kb, 5-150 Mb. Released on/after 2025-05-01 preferred; older non-BFD genomes fill up quotas.
Adds columns species_in_bfd / genus_in_bfd so the benchmark can tell "known species" from "novel species".
"""
import sys, random, json, collections

d, quota, out = sys.argv[1], json.loads(sys.argv[2]), sys.argv[3]
bfd = [l.split("\t") for l in open(f"{d}/bfd_meta.tsv").read().splitlines()[1:]]
bfd_nums = {r[0].split("_")[1].split(".")[0] for r in bfd}
bfd_species = {r[7] for r in bfd if r[7]}
bfd_genera = {r[6] for r in bfd if r[6]}

lin = {}
with open(f"{d}/rankedlineage.dmp") as fh:
    for line in fh:
        f = [x.strip() for x in line.rstrip("\n").split("\t|")]
        # tax_id name species genus family order class phylum kingdom superkingdom
        lin[f[0]] = dict(species=f[2], genus=f[3], family=f[4], order=f[5], cls=f[6], phylum=f[7])

rows = [l.split("\t") for l in open(f"{d}/fungi_all.tsv").read().splitlines()[1:]]
cand = []
for r in rows:
    acc, name, taxid, rel, level, tlen, n50 = r[:7]
    if acc.split("_")[1].split(".")[0] in bfd_nums:
        continue
    try:
        tlen, n50 = int(tlen), int(n50)
    except ValueError:
        continue
    if n50 < 20000 or not (5e6 <= tlen <= 150e6):
        continue
    l = lin.get(taxid)
    if not l or not l["phylum"]:
        continue
    cand.append(dict(acc=acc, name=name, taxid=taxid, rel=rel, level=level, tlen=tlen, n50=n50, **l))

random.seed(1)
chosen = []
by_ph = collections.defaultdict(list)
for c in cand:
    by_ph[c["phylum"]].append(c)
for ph, items in by_ph.items():
    q = quota.get(ph, quota.get("other", 0))
    # prefer recent releases, then random, and spread over species
    random.shuffle(items)
    items.sort(key=lambda c: c["rel"] < "2025-05-01")
    seen_sp = collections.Counter()
    pick = []
    for c in items:
        if len(pick) >= q:
            break
        if seen_sp[c["species"]] >= 4:
            continue
        seen_sp[c["species"]] += 1
        pick.append(c)
    chosen += pick
with open(out + ".meta.tsv", "w") as m:
    m.write("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\tspecies_in_bfd\tgenus_in_bfd\trelease\tn50\ttotal_len\n")
    for c in chosen:
        m.write("\t".join([c["acc"], c["taxid"], c["phylum"], c["cls"], c["order"], c["family"], c["genus"], c["species"],
                           str(c["species"] in bfd_species), str(c["genus"] in bfd_genera), c["rel"], str(c["n50"]), str(c["tlen"])]) + "\n")
open(out + ".accessions.txt", "w").write("\n".join(c["acc"] for c in chosen) + "\n")
cnt = collections.Counter((c["phylum"], c["species"] in bfd_species) for c in chosen)
print("candidates", len(cand), "chosen", len(chosen))
for (ph, known), n in sorted(cnt.items()):
    print(f"  {ph:22s} species_in_bfd={known!s:5s} {n}")
print("genus_in_bfd False:", sum(1 for c in chosen if c["genus"] not in bfd_genera))
