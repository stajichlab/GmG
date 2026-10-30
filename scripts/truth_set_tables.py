#!/usr/bin/env python3
"""Tables for reports/validation/truth_set.md (#39) from the validate_run.py outputs and the disagreement drafts.

usage: truth_set_tables.py SCORE_DIR OUT_TABLES.md DISAGREEMENTS_OUT.tsv [PREVIOUS_DISAGREEMENTS.tsv]
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
V, VN = "known.i3", "novel.i3"  # production index 2026.10.3, default and novel thresholds
M3 = [("loo", "leave-one-out"), ("los", "leave-species-out"), ("log", "leave-genus-out"), ("indep", "outside the index (21)")]
md.append("### T1. All inputs, called fraction / precision [95% Wilson] (calls made), index 2026.10.3\n\n" + table(M3, variant=V, kind="raw"))
for grp, mode, title in [("species_in_bfd=True", "loo", "leave-one-out, species still in the reference (answerable)"),
                         ("species_in_bfd=False", "loo", "leave-one-out, species a singleton in the reference (not answerable)"),
                         ("evidence=type_material", "loo", "leave-one-out, type material only"),
                         ("evidence=reference_genome", "loo", "leave-one-out, reference genomes (not type)"),
                         ("genus_in_bfd=True", "los", "leave-species-out, genus still in the reference"),
                         ("genus_in_bfd=False", "los", "leave-species-out, genus is a singleton (also absent)"),
                         ("genus_in_bfd=False", "log", "leave-genus-out, genus absent")]:
    t = table([(mode, title)], group=grp, variant=V, kind="raw")
    if t.count("\n") > 2:
        md.append(f"### {title}\n\n" + t)
md.append("### T2. `--thresholds novel`, index 2026.10.3\n\n" + table([("los", "leave-species-out"), ("log", "leave-genus-out"), ("indep", "outside the index")], variant=VN, kind="raw"))
M3b = M3[:3]
before = [(m, l + ", 2026.10.2 (43 audit genomes kept)") for m, l in M3b]
after = [(m, l + ", 2026.10.3 (57 left out)") for m, l in M3b]
hid = [(m, l + ", 2026.10.2 with the 43 hidden") for m, l in M3b]
t_b = table(before, variant="known", kind="raw").splitlines()
t_h = table(hid, variant="known.noaudit", kind="raw").splitlines()
t_a = table(after, variant=V, kind="raw").splitlines()
rows = [t_b[2 + i] for i in range(3)]
md.append("### T4. Before and after leaving the 57 ambiguous-label genomes out\n\n" + "\n".join(t_b[:2] + [x for i in range(3) for x in (t_b[2 + i], t_h[2 + i], t_a[2 + i])]))
nod = [(m, f"{l}, 2026.10.1 (no decoys)") for m, l in M3b]
dec_raw = [(m, f"{l}, 2026.10.2 (decoys, raw)") for m, l in M3b]
dec_aw = [(m, f"{l}, 2026.10.2 (decoys, flag-aware)") for m, l in M3b]
t1 = table(nod, variant="known.nodecoy1", kind="raw").splitlines()
t2 = table(dec_raw, variant="known", kind="raw").splitlines()
t3 = table(dec_aw, variant="known", kind="aware").splitlines()
t4 = table([(m, f"{l}, 2026.10.3 (decoys, call re-made without them)") for m, l in M3b], variant=V, kind="raw").splitlines()
md.append("### T3. Decoys: phylum to genus, same taxonomy\n\n" + "\n".join(t1[:2] + [x for i in range(3) for x in (t1[2 + i], t2[2 + i], t3[2 + i], t4[2 + i])]))

# disagreements: combined log and category tables (production run; dataset names without the index suffix)
allrows = []
for ds in ("loo.known", "los.known", "log.known", "indep.known", "los.novel", "log.novel", "indep.novel"):
    p = os.path.join(sd, f"disagreements_draft.{ds}.i3.tsv")
    if not os.path.exists(p):
        continue
    for r in csv.DictReader(open(p), delimiter="\t"):
        r["dataset"] = ds
        allrows.append(r)
prev = {}
if len(sys.argv) > 4 and os.path.exists(sys.argv[4]):
    for r in csv.DictReader(open(sys.argv[4]), delimiter="\t"):
        prev[(r["dataset"], r["id"], r["rank"], r["call"])] = r["category"]
for r in allrows:
    r["prev_category"] = prev.get((r["dataset"], r["id"], r["rank"], r["call"]), "(not a disagreement in 2026.10.2)")
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
# which rows changed category since the previous log
ch = collections.Counter()
for r in allrows:
    if r["dataset"].endswith(".known"):
        ch[(r["prev_category"] == r["category"], r["prev_category"] if r["prev_category"] != r["category"] else "")] += 1
tr = collections.Counter((r["prev_category"], r["category"]) for r in allrows if r["dataset"].endswith(".known") and r["prev_category"] != r["category"])
now = {(r["dataset"], r["id"], r["rank"], r["call"]) for r in allrows}
gone = collections.Counter(c for k, c in prev.items() if k[0].endswith(".known") and k not in now)
same = sum(1 for r in allrows if r["dataset"].endswith(".known") and r["prev_category"] == r["category"])
both = sum(1 for r in allrows if r["dataset"].endswith(".known") and not r["prev_category"].startswith("(not"))
md.append("### T6. Rows whose category changed against the 2026.10.2 log (known thresholds)\n\n"
          f"Rows in both logs: {both} ({same} with the same category). Rows only in the 2026.10.3 log: "
          f"{sum(1 for r in allrows if r['dataset'].endswith('.known') and r['prev_category'].startswith('(not'))}. "
          f"Rows only in the 2026.10.2 log (no longer a disagreement): {sum(gone.values())} "
          f"({', '.join(f'{k} {v}' for k, v in gone.most_common())}).\n\n| 2026.10.2 | 2026.10.3 | rows |\n|---|---|---|\n" +
          "".join(f"| {a} | {b} | {n} |\n" for (a, b), n in sorted(tr.items(), key=lambda kv: -kv[1])))
open(out_md, "w").write("\n".join(md))

cols = ["dataset", "id", "rank", "truth", "call", "category", "prev_category", "evidence", "decision", "date", "status", "ani", "top_hit", "top_hit_species", "bfd_label", "audit", "n50", "qc_flags", "reviewed"]
with open(out_dis, "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in allrows:
        if r["decision"] == "draft (auto)":
            r["decision"] = "rule-based (not reviewed individually)"
        o.write("\t".join(r[c].replace("\t", " ") for c in cols) + "\n")
print(len(allrows), "disagreement rows written")
