#!/usr/bin/env python3
"""Trusted truth set (#39): NCBI RefSeq fungal assemblies, with type-material / reference flags, their NCBI lineage and the
reference genome each one corresponds to in the index.

usage: make_truth_set.py REFSEQ.jsonl INDEX_META.tsv LINEAGE.tsv OUT_PREFIX
  REFSEQ.jsonl   `datasets summary genome taxon fungi --assembly-source RefSeq --as-json-lines`
  INDEX_META.tsv release meta_full.tsv (index ids; the numeric part of the accession decides "same assembly")
  LINEAGE.tsv    normalize_taxonomy.py output for the truth taxids (id = RefSeq accession): pinned NCBI lineage = the truth
Writes OUT_PREFIX.all.tsv: every assembly with its truth lineage, evidence, level, n50, total_len, bfd_id (the index genome that
is the same assembly, if any) and how many other genomes of its species / genus the index holds once that genome is left out.
"""
import csv, json, sys, collections

jl, meta, lin, out = sys.argv[1:5]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def num(acc):
    return acc.split("_")[1].split(".")[0]


idx = list(csv.DictReader(open(meta), delimiter="\t"))
by_num = collections.defaultdict(list)
for r in idx:
    by_num[num(r["id"])].append(r["id"])
sp_n = collections.Counter(r["species"] for r in idx if r["species"])
ge_n = collections.Counter(r["genus"] for r in idx if r["genus"])
truth = {r["id"]: r for r in csv.DictReader(open(lin), delimiter="\t")}

rows = []
for line in open(jl):
    j = json.loads(line)
    acc = j["accession"]
    ai, st, org = j.get("assembly_info", {}), j.get("assembly_stats", {}), j.get("organism", {})
    tm = j.get("type_material") or {}
    cat = (ai.get("refseq_category") or "").lower()
    ev = "type_material" if tm else ("reference_genome" if cat in ("reference genome", "representative genome") else "refseq")
    t = truth.get(acc)
    if not t or not t["phylum"]:
        continue
    paired = j.get("paired_accession") or ai.get("paired_assembly", {}).get("accession") or ""
    cand = by_num.get(num(acc), []) + (by_num.get(num(paired), []) if paired else [])
    cand = list(dict.fromkeys(cand))
    bfd_id = cand[0] if cand else ""
    sp, ge = t["species"], t["genus"]
    n_sp = (sp_n[sp] - (1 if bfd_id else 0)) if sp else 0
    n_ge = (ge_n[ge] - (1 if bfd_id else 0)) if ge else 0
    rows.append(dict(id=acc, taxid=t["taxid"], **{r: t[r] for r in RANKS}, organism=org.get("organism_name", ""),
                     type_material=("yes" if tm else "no"), type_label=tm.get("type_label", ""), refseq_category=cat,
                     evidence=ev, level=ai.get("assembly_level", ""), n50=st.get("contig_n50", ""),
                     total_len=st.get("total_sequence_length", ""), release=ai.get("release_date", ""),
                     bfd_id=bfd_id, bfd_ids_all=",".join(cand), in_bfd=("yes" if bfd_id else "no"), species_others_in_ref=n_sp, genus_others_in_ref=n_ge))
cols = list(rows[0])
with open(out + ".all.tsv", "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: r["id"]):
        o.write("\t".join(str(r[c]) for c in cols) + "\n")
print(len(rows), "assemblies;", sum(r["in_bfd"] == "yes" for r in rows), "in the index;",
      sum(r["type_material"] == "yes" for r in rows), "type material;",
      sum(r["type_material"] == "yes" and r["in_bfd"] == "yes" for r in rows), "type material in the index")
