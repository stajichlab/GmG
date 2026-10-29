#!/usr/bin/env python3
"""Score a `gmg run` results.tsv against provider names ("Genus species ...", strain labels) via the NCBI taxdump (#40).

usage: score_by_name.py RESULTS.tsv LABELS.tsv TAXDUMP_DIR OUT_PREFIX     (LABELS: columns sample, name)
Truth lineage = the NCBI lineage of the name's genus (first word; brackets and "cf." removed), found in names.dmp as a
scientific name or synonym, so old names (Hypocrea, Cochliobolus) resolve to the current genus. Names whose genus is not
in the taxdump, or that are not genus-first ("Uncultured ..."), are listed as unscored. A call is correct at a rank when
it equals the truth lineage name at that rank. Writes OUT_PREFIX.tsv (one row per input; keep it private if the names are)
and OUT_PREFIX.md (aggregates only). Species is not scored: provider species names are not curated.
"""
import sys, re, csv, collections
res_p, lab_p, td, out = sys.argv[1:5]
RANKS = ["phylum", "class", "order", "family", "genus"]
labels = {r["sample"]: r["name"] for r in csv.DictReader(open(lab_p), delimiter="\t")}
want = set()
genus_of = {}
for s, n in labels.items():
    n = re.sub(r"[\[\]]|\bcf\.?\b|\baff\.?\b", "", n).strip()
    w = n.split()[0].split("/")[0] if n else ""
    if w and w[0].isupper(): genus_of[s] = w; want.add(w)
by_name = collections.defaultdict(set)   # name -> taxids
for l in open(f"{td}/names.dmp"):
    f = [x.strip() for x in l.split("|")]
    if f[1] in want and f[3] in ("scientific name", "synonym", "includes", "equivalent name"): by_name[f[1]].add(int(f[0]))
parent, rank = {}, {}
for l in open(f"{td}/nodes.dmp"):
    f = [x.strip() for x in l.split("|")]
    parent[int(f[0])] = int(f[1]); rank[int(f[0])] = f[2]
need = {t for ts in by_name.values() for t in ts}
anc = {}
for t in need:
    x, d = t, {}
    while x in parent and parent[x] != x:
        d.setdefault(rank[x], x); x = parent[x]
    anc[t] = d
sci = {}
allt = {v for d in anc.values() for v in d.values()}
for l in open(f"{td}/names.dmp"):
    f = [x.strip() for x in l.split("|")]
    if f[3] == "scientific name" and int(f[0]) in allt: sci[int(f[0])] = f[1]
def truth(g):
    best = None
    for t in by_name.get(g, ()):
        d = anc[t]
        if rank[t] == "genus" or "genus" in d:
            lin = {r: sci.get(d.get(r), "") for r in RANKS}
            if best is None or len(lin["phylum"]) > len(best["phylum"]): best = lin
    return best
rows = list(csv.DictReader(open(res_p), delimiter="\t"))
tab = []; stat = collections.Counter(); unscored = 0
for r in rows:
    s = r["query_id"]; stat[r["status"]] += 1
    g = genus_of.get(s); tr = truth(g) if g else None
    if not tr: unscored += 1; tab.append(dict(sample=s, provider_genus=g or "", status=r["status"])); continue
    d = dict(sample=s, provider_genus=g, status=r["status"])
    for k in RANKS:
        c = r.get(f"call_{k}", ""); d[f"truth_{k}"] = tr[k]; d[f"call_{k}"] = c
        d[f"ok_{k}"] = "" if not c or not tr[k] else str(c == tr[k])
    tab.append(d)
cols = ["sample", "provider_genus", "status"] + [f"{p}_{k}" for k in RANKS for p in ("truth", "call", "ok")]
w = csv.DictWriter(open(out + ".tsv", "w"), fieldnames=cols, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(tab)
md = [f"# Scored by provider genus name ({len(rows)} inputs)\n", f"status: {dict(stat)}; {unscored} without a usable genus name in the taxdump.\n",
      "| rank | with truth | called | correct | precision | called share |", "|---|---|---|---|---|---|"]
scored_rows = [d for d in tab if "truth_phylum" in d]
for k in RANKS:
    wt = [d for d in scored_rows if d[f"truth_{k}"]]
    called = [d for d in wt if d[f"call_{k}"]]
    ok = sum(d[f"ok_{k}"] == "True" for d in called)
    md.append(f"| {k} | {len(wt)} | {len(called)} | {ok} | {ok/len(called):.3f} | {len(called)/len(wt):.3f} |" if called else f"| {k} | {len(wt)} | 0 | 0 | - | 0 |")
wrong = collections.Counter((d["truth_phylum"], d["call_phylum"]) for d in scored_rows if d["ok_phylum"] == "False")
md.append("\nWrong phylum calls (truth, call): " + (", ".join(f"{a}->{b} x{n}" for (a, b), n in wrong.most_common(8)) or "none"))
gw = [d for d in scored_rows if d["ok_genus"] == "False"]
md.append(f"\nWrong genus calls: {len(gw)}. Most frequent (provider genus -> call): " + ", ".join(f"{d['provider_genus']}->{d['call_genus']}" for d in gw[:15]))
open(out + ".md", "w").write("\n".join(md) + "\n"); print("\n".join(md))
