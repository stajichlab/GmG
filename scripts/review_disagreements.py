#!/usr/bin/env python3
"""List every disagreement between a `gmg run` results.tsv and a truth table, with the evidence needed to review it (#39).

usage: review_disagreements.py DATASET RESULTS.tsv[,RESULTS.tsv...] TRUTH.tsv TRUTH_BFDLABEL.tsv|- INDEX_META.tsv TAXDUMP_DIR AUDIT_DIR OUT.tsv
                               [--reviews REVIEWS.tsv] [--excluded-rank species|genus]
A disagreement is a CALL that differs from the truth at a rank where the truth is known; a withheld rank is not one.
Draft categories are rule-based (the table says which rule fired); REVIEWS.tsv (dataset, id, rank, category, decision) holds the
manual verdicts and overrides the draft. Categories: tool_error, reference_label_error, truth_label_error, species_complex,
taxonomy_synonym, contamination, low_quality_input. Columns of OUT: dataset id rank truth call category evidence decision date
plus the numbers behind them (status, ani, top hit, BFD label of the query, audit lists, quality).
--excluded-rank: the reference holds nothing of this rank for the query (leave-species/genus-out), so a call there cannot be right.
"""
import csv, collections, datetime, os, sys

args = [x for x in sys.argv[1:]]
opt = {}
for k in ("--reviews", "--excluded-rank"):
    if k in args:
        i = args.index(k); opt[k] = args[i + 1]; del args[i:i + 2]
dataset, res_p, truth_p, bfd_p, meta_p, taxdump, audit, out = args
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
res = [r for p in res_p.split(",") for r in rd(p)]
truth = {r["id"]: r for r in rd(truth_p)}
bfd = {r["id"]: r for r in rd(bfd_p)} if bfd_p != "-" else {}
meta = {r["id"]: r for r in rd(meta_p)}
names = collections.defaultdict(set)  # any NCBI name -> taxids (synonyms resolve to the same taxid)
with open(os.path.join(taxdump, "names.dmp")) as f:
    for l in f:
        p = l.rstrip("\n").split("\t|\t")
        names[p[1].lower()].add(p[0])
audit_sets = {}
for name in ("exclude_ids.txt", "loo_suspects.tsv", "hub_references.tsv", "size_outliers.tsv", "log_conflicts.tsv"):
    pth = os.path.join(audit, name)
    if os.path.exists(pth):
        audit_sets[name.split(".")[0]] = {l.split("\t")[0].strip() for l in open(pth) if l.strip() and not l.startswith(("id", "top_id", "#"))}
reviews = {}
if "--reviews" in opt and os.path.exists(opt["--reviews"]):
    for r in rd(opt["--reviews"]):
        reviews[(r["dataset"], r["id"], r["rank"])] = r


def synonym(a, b):
    return bool(a and b and names.get(a.lower()) and names.get(a.lower()) & names.get(b.lower(), set()))


def audit_hits(i):
    return ",".join(k for k, v in audit_sets.items() if i in v)


rows = []
for r in res:
    t = truth.get(r["query_id"])
    if not t or r["status"] == "error":
        continue
    own = bfd.get(r["query_id"], {})
    top = r.get("top_hit", "")
    for rank in RANKS:
        tv, cv = t[rank], r[f"call_{rank}"]
        if not tv or not cv or tv == cv:
            continue
        flags = r.get("qc_flags", "")
        n50 = int(r["n50"]) if r.get("n50", "").isdigit() else 0
        nh = int(r["n_hashes"]) if r.get("n_hashes", "").isdigit() else 0
        topmeta = meta.get(top, {})
        why = []
        if synonym(tv, cv):
            cat, why = "taxonomy_synonym", [f"'{tv}' and '{cv}' share an NCBI taxid"]
        elif n50 < 10000 or nh < 1000 or any(f in flags for f in ("fragmented", "short_genome", "few_hashes")):
            cat, why = "low_quality_input", [f"n50 {n50}, {nh} hashes, flags {flags or '-'}"]
        elif rank == "species" and t["genus"] == r["call_genus"] and r["call_genus"]:
            cat, why = "species_complex", [f"same genus {t['genus']}; ANI {r.get('ani', '')}; confusable congeners"]
        elif audit_hits(top) and cv == topmeta.get(rank, ""):
            cat, why = "reference_label_error", [f"top hit {top} is on the label-audit list ({audit_hits(top)})"]
        elif own and own.get(rank) == cv and own.get(rank) != tv:
            cat, why = "truth_label_error", [f"the index label of this genome ({own.get(rank)}) agrees with the call, NCBI RefSeq says {tv}"]
        else:
            cat = "tool_error"
            if opt.get("--excluded-rank") == rank:
                why = [f"call at the excluded rank ({rank}) the reference cannot hold (overcall)"]
        rv = reviews.get((dataset, r["query_id"], rank))
        rows.append(dict(
            dataset=dataset, id=r["query_id"], rank=rank, truth=tv, call=cv,
            category=rv["category"] if rv else cat, auto_category=cat,
            evidence="; ".join(why), decision=rv["decision"] if rv else "draft (auto)",
            date=datetime.date.today().isoformat(), status=r["status"], ani=r.get("ani", ""), top_hit=top,
            top_hit_species=r.get("top_hit_species", ""), bfd_label=own.get(rank, ""), audit=audit_hits(top) or audit_hits(r["query_id"]),
            n50=n50, n_hashes=nh, qc_flags=flags, reviewed="yes" if rv else "no"))
cols = list(rows[0]) if rows else ["dataset", "id", "rank"]
with open(out, "w") as o:
    o.write("\t".join(cols) + "\n")
    for x in rows:
        o.write("\t".join(str(x[c]) for c in cols) + "\n")
c = collections.Counter(x["category"] for x in rows)
print(dataset, len(rows), "disagreements:", dict(c))
