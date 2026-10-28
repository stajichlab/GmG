#!/usr/bin/env python3
"""Pick per-species representatives with the BFD rule: highest BUSCO genome completeness, then highest N50, then ASMID.

usage: pick_reps.py META.tsv BFD_TABLES_DIR OUT_PREFIX [K ...]
writes OUT_PREFIX.order.txt (all ids, best first), OUT_PREFIX.reps<K>.txt (ids) and OUT_PREFIX.meta_reps<K>.tsv (meta rows) for each K (top K per species).
Genomes with an empty species are each kept (they cannot be grouped).
"""
import sys, duckdb, csv

meta, tdir, out = sys.argv[1:4]
ks = [int(k) for k in sys.argv[4:]] or [1, 3]
c = duckdb.connect()
q = c.execute(f"""
  select a.ASMID, coalesce(b.complete_pct, -1) as busco, coalesce(a.N50_bp, 0) as n50
  from read_parquet('{tdir}/asm_stats.parquet') a left join read_parquet('{tdir}/busco_genome.parquet') b using (ASMID)
""").fetchall()
qual = {r[0]: (r[1], r[2]) for r in q}
rows = list(csv.reader(open(meta), delimiter="\t"))
header, rows = rows[0], rows[1:]
# priority order of all genomes (BFD rule) for `gmg cluster --order`
allrows = sorted(rows, key=lambda r: (-qual.get(r[0], (-1, 0))[0], -qual.get(r[0], (-1, 0))[1], r[0]))
open(f"{out}.order.txt", "w").write("\n".join(r[0] for r in allrows) + "\n")
by_sp = {}
for r in rows:
    by_sp.setdefault(r[7] or f"_nosp_{r[0]}", []).append(r)
for k in ks:
    keep = []
    for sp, g in by_sp.items():
        g.sort(key=lambda r: (-qual.get(r[0], (-1, 0))[0], -qual.get(r[0], (-1, 0))[1], r[0]))
        keep += g[:k]
    keep.sort(key=lambda r: r[0])
    open(f"{out}.reps{k}.txt", "w").write("\n".join(r[0] for r in keep) + "\n")
    with open(f"{out}.meta_reps{k}.tsv", "w") as f:
        f.write("\t".join(header) + "\n")
        for r in keep:
            f.write("\t".join(r) + "\n")
    print(f"K={k}: {len(keep)} reference genomes for {len(by_sp)} species (no quality row: {sum(1 for r in keep if r[0] not in qual)})")
