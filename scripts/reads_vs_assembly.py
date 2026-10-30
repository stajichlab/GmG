#!/usr/bin/env python3
"""Agreement of the read-based and the assembly-based call for the same strain (#42), independent of the lab's label.

usage: reads_vs_assembly.py READS_RESULTS.tsv ASSEMBLY_RESULTS.tsv TRUTH.tsv
READS_RESULTS: `gmg run` of read sets whose ids are <Genus_species_STRAIN>; ASSEMBLY_RESULTS: `gmg run` of assemblies whose ids are the strain
(e.g. NRRL_1807); TRUTH.tsv: the lab-label lineage of the reads (name_truth.py). A read set is matched to the assembly whose id is the end of the read id.
Prints: matched pairs; agreement of the two calls at phylum and genus (both called); and, where the read call disagrees with the lab label, how often the
assembly makes the same call as the reads (then the label, not the reads method, is the odd one out).
"""
import csv, sys

reads = {r["query_id"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
asm = {r["query_id"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
tru = {r["id"]: r for r in csv.DictReader(open(sys.argv[3]), delimiter="\t")}
pairs = []
for rid, r in reads.items():
    cands = [a for a in asm if rid.endswith("_" + a)]
    if len(cands) >= 1:
        pairs.append((r, asm[max(cands, key=len)], tru.get(rid)))
print(f"{len(pairs)} of {len(reads)} read sets matched to an assembly of the same strain")
for rank in ("phylum", "genus"):
    both = [(r, a) for r, a, t in pairs if r[f"call_{rank}"] and a[f"call_{rank}"]]
    agree = sum(r[f"call_{rank}"] == a[f"call_{rank}"] for r, a in both)
    print(f"{rank}: both called {len(both)}, same call {agree} ({agree / max(1, len(both)):.3f})")
for rank in ("phylum", "genus"):
    dis = [(r, a, t) for r, a, t in pairs if t and r[f"call_{rank}"] and r[f"call_{rank}"] != t[rank]]
    same = [(r, a, t) for r, a, t in dis if a[f"call_{rank}"] == r[f"call_{rank}"]]
    print(f"{rank}: read call differs from the lab label for {len(dis)}; the assembly makes the same call as the reads for {len(same)}")
    if rank == "genus":
        for r, a, t in dis[:40]:
            print("   ", t["label_genus"][:18], "| reads:", r["call_genus"][:18], r["status"][:4], "| assembly:", (a["call_genus"] or "-")[:18], a["status"][:4])
