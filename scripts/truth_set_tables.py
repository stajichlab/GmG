#!/usr/bin/env python3
"""Tables for reports/validation/truth_set.md (#39) from the validate_run.py outputs and the disagreement drafts.

usage: truth_set_tables.py SCORE_DIR OUT_TABLES.md DISAGREEMENTS_OUT.tsv
SCORE_DIR has MODE.VARIANT.{raw,aware}.refseq.tsv (validate_run.py), disagreements_draft.MODE.VARIANT.tsv (review_disagreements.py).
Writes the Markdown tables and the combined disagreement log (data/validation/disagreements.tsv).
"""
import csv, math, os, sys, collections

sd, out_md, out_dis = sys.argv[1:4]
RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def load(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    d = collections.defaultdict(dict)
    for r in rows:
        d[r["group"]][r["rank"]] = r
    return d


def cell(r):
    if r is None:
        return "n/a"
    n, called = int(r["n"]), float(r["called"])
    nc = round(called * n)
    prec = r["precision"]
    if nc == 0 or prec == "":
        return f"{called:.2f} / n/a ({nc})"
    k = round(float(prec) * nc)
    lo, hi = wilson(k, nc)
    return f"{called:.2f} / {float(prec):.3f} [{lo:.3f}, {hi:.3f}] ({nc})"


def table(modes, group="all", variant="known", kind="aware"):
    head = "| mode | n | " + " | ".join(RANKS) + " |\n|---|---|" + "---|" * len(RANKS) + "\n"
    body = ""
    for m, label in modes:
        p = os.path.join(sd, f"{m}.{variant}.{kind}.refseq.tsv")
        if not os.path.exists(p):
            continue
        g = load(p).get(group)
        if not g:
            continue
        n = g["phylum"]["n"]
        body += f"| {label} | {n} | " + " | ".join(cell(g.get(r)) for r in RANKS) + " |\n"
    return head + body


md = []
M3 = [("loo", "leave-one-out"), ("los", "leave-species-out"), ("log", "leave-genus-out"), ("indep", "outside the index (21)")]
md.append("### T1. All inputs, called fraction / precision [95% Wilson] (calls made), flag-aware\n\n" + table(M3))
md.append("### T1b. Same, raw (a distant decoy phylum call counted as a call)\n\n" + table(M3, kind="raw"))
for grp, mode, title in [("species_in_bfd=True", "loo", "leave-one-out, species still in the reference (answerable)"),
                         ("species_in_bfd=False", "loo", "leave-one-out, species a singleton in the reference (not answerable)"),
                         ("evidence=type_material", "loo", "leave-one-out, type material only"),
                         ("evidence=reference_genome", "loo", "leave-one-out, reference genomes (not type)"),
                         ("genus_in_bfd=True", "los", "leave-species-out, genus still in the reference"),
                         ("genus_in_bfd=False", "los", "leave-species-out, genus is a singleton (also absent)"),
                         ("genus_in_bfd=True", "log", "leave-genus-out, genus still present (should be empty)"),
                         ("genus_in_bfd=False", "log", "leave-genus-out, genus absent")]:
    t = table([(mode, title)], group=grp)
    if t.count("\n") > 2:
        md.append(f"### {title}\n\n" + t)
md.append("### T2. `--thresholds novel`, flag-aware\n\n" + table([("los", "leave-species-out"), ("log", "leave-genus-out"), ("indep", "outside the index")], variant="novel"))
nod = [(m, f"{l} without decoys") for m, l in M3[:3]]
md.append("### T3. Index 2026.10.1 (same taxonomy, no decoys), raw\n\n" + table(nod, variant="known.nodecoy1", kind="raw") +
          "\n### T3b. Index 2026.10.2 (with decoys), raw\n\n" + table(M3[:3], kind="raw"))
md.append("### T4. Reference genomes on the label-audit list (43) hidden, flag-aware\n\n" + table([(m, l + ", audit list hidden") for m, l in M3[:3]], variant="known.noaudit") +
          "\n### T4b. Audit list kept (the production index), flag-aware\n\n" + table(M3[:3]))

# disagreements: combined log and category tables
allrows, tabs = [], []
for ds in ("loo.known", "los.known", "log.known", "indep.known", "los.novel", "log.novel", "indep.novel"):
    p = os.path.join(sd, f"disagreements_draft.{ds}.tsv")
    if not os.path.exists(p):
        continue
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    allrows += rows
cats = ["tool_error", "reference_label_error", "taxonomy_synonym", "species_complex", "low_quality_input", "truth_label_error", "contamination"]
head = "| dataset | disagreements | " + " | ".join(cats) + " | status ok (confident) | reviewed by hand |\n|---|---|" + "---|" * (len(cats) + 2) + "\n"
body = ""
for ds in ("loo.known", "los.known", "log.known", "los.novel", "log.novel", "indep.known"):
    rr = [r for r in allrows if r["dataset"] == ds]
    c = collections.Counter(r["category"] for r in rr)
    ok = sum(r["status"] == "ok" for r in rr)
    rv = sum(r["reviewed"] == "yes" for r in rr)
    body += f"| {ds} | {len(rr)} | " + " | ".join(str(c.get(k, 0)) for k in cats) + f" | {ok} | {rv} |\n"
md.append("### T5. Disagreements (a wrong call at a rank with a known truth) by category\n\n" + head + body)
# by rank x category for the three modes
for ds in ("loo.known", "los.known", "log.known"):
    rr = [r for r in allrows if r["dataset"] == ds]
    c = collections.Counter((r["rank"], r["category"]) for r in rr)
    head = "| rank | " + " | ".join(cats[:5]) + " | total | of which status ok |\n|---|" + "---|" * 7 + "\n"
    body = ""
    for rk in RANKS:
        tot = sum(c.get((rk, k), 0) for k in cats)
        ok = sum(1 for r in rr if r["rank"] == rk and r["status"] == "ok")
        body += f"| {rk} | " + " | ".join(str(c.get((rk, k), 0)) for k in cats[:5]) + f" | {tot} | {ok} |\n"
    md.append(f"### {ds}: by rank\n\n" + head + body)
open(out_md, "w").write("\n".join(md))

cols = ["dataset", "id", "rank", "truth", "call", "category", "evidence", "decision", "date", "status", "ani", "top_hit", "top_hit_species", "bfd_label", "audit", "n50", "qc_flags", "reviewed"]
with open(out_dis, "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in allrows:
        if r["decision"] == "draft (auto)":
            r["decision"] = "rule-based (not reviewed individually)"
        o.write("\t".join(r[c].replace("\t", " ") for c in cols) + "\n")
print(len(allrows), "disagreement rows written")
