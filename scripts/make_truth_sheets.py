#!/usr/bin/env python3
"""Sample sheets and truth tables for the trusted truth set (#39).

usage: make_truth_sheets.py TRUTH.all.tsv INDEX_META.tsv OUT_DIR [--shard N]
For every assembly in the truth set that is in the index (bfd_id) the sheets for leave-one-out (loo), leave-species-out (los) and
leave-genus-out (log) runs, in shards of N rows (default 50): columns id, path (../genomes/ACC.fna.gz), expected (truth species or
genus text, used when the lineage of `ref` is not in the index), ref (the index genome that is this assembly). The assemblies that
are not in the index form `indep` (no exclusion). Truth tables per mode (truth_MODE.tsv, NCBI lineage of the RefSeq accession) and
`truth_bfdlabel_MODE.tsv` (the lineage the index holds for the same genome: "BFD label" truth) carry species_in_bfd / genus_in_bfd /
family_in_bfd as the query sees the reference after the exclusion, so validate_run.py can report the false-call rates.
"""
import csv, collections, os, sys

truth_p, meta_p, out = sys.argv[1:4]
shard = int(sys.argv[sys.argv.index("--shard") + 1]) if "--shard" in sys.argv else 50
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
os.makedirs(out, exist_ok=True)
rows = list(csv.DictReader(open(truth_p), delimiter="\t"))
meta = {r["id"]: r for r in csv.DictReader(open(meta_p), delimiter="\t")}
ids_by = {r: collections.defaultdict(set) for r in ("species", "genus", "family")}
for m in meta.values():
    for r in ids_by:
        if m[r]:
            ids_by[r][m[r]].add(m["id"])
EX = {"loo": "self", "los": "species", "log": "genus"}


def hidden_for(mode, own):
    """index genomes the query cannot see: itself, its species, or its genus (as `gmg run --exclude-level` hides them)"""
    if own is None:
        return set()
    h = {own["id"]}
    if mode == "los" and own["species"]:
        h |= ids_by["species"][own["species"]]
    if mode == "log" and own["genus"]:
        h |= ids_by["genus"][own["genus"]]
    return h


def remaining(rank, taxon, hidden):
    """number of index genomes of `taxon` (at `rank`) the query still sees"""
    return len(ids_by[rank][taxon] - hidden) if taxon else 0


def write_truth(path, recs, lineage_of, mode):
    cols = ["id", *RANKS, "species_in_bfd", "genus_in_bfd", "family_in_bfd", "type_material", "evidence", "level", "n50", "total_len"]
    with open(path, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in recs:
            lin = lineage_of(r)
            own = meta.get(r["bfd_id"]) if r["bfd_id"] else None
            hid = hidden_for(mode, own)
            vals = [r["id"], *[lin[x] for x in RANKS],
                    str(remaining("species", lin["species"], hid) > 0), str(remaining("genus", lin["genus"], hid) > 0),
                    str(remaining("family", lin["family"], hid) > 0), r["type_material"], r["evidence"], r["level"], r["n50"], r["total_len"]]
            o.write("\t".join(vals) + "\n")


inb = [r for r in rows if r["bfd_id"]]
outb = [r for r in rows if not r["bfd_id"]]
for mode in ("loo", "los", "log"):
    recs = [r for r in inb]
    for i in range(0, len(recs), shard):
        with open(f"{out}/{mode}.{i // shard:03d}.tsv", "w") as o:
            o.write("id\tpath\texpected\tref\n")
            for r in recs[i:i + shard]:
                exp = r["species"] if mode != "log" else (r["genus"] or r["species"])
                o.write(f"{r['id']}\t../genomes/{r['id']}.fna.gz\t{exp}\t{r['bfd_id']}\n")
    write_truth(f"{out}/truth_{mode}.tsv", recs, lambda r: r, mode)
    write_truth(f"{out}/truth_bfdlabel_{mode}.tsv", recs, lambda r: meta[r["bfd_id"]], mode)
with open(f"{out}/indep.000.tsv", "w") as o:
    o.write("id\tpath\texpected\tref\n")
    for r in outb:
        o.write(f"{r['id']}\t../genomes/{r['id']}.fna.gz\t{r['species']}\t\n")
write_truth(f"{out}/truth_indep.tsv", outb, lambda r: r, "none")
print(len(inb), "in the index,", len(outb), "outside;", (len(inb) + shard - 1) // shard, "shards per mode")
