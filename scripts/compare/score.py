#!/usr/bin/env python3
"""Score the tools of the comparison (#46) with one code path. usage: score.py COMPARE_DIR REPO_DIR OUT.tsv [META.tsv]
For every tool and query set it turns the tool's output into a ranked hit list (id, ANI), drops the hits the query may not see (itself; its
species; its genus: leave-one/species/genus-out, applied after the search, which is valid for nearest-genome tools), and makes calls from the
best remaining hit in two ways:
  raw     the lineage of the best hit at every rank (what "nearest genome" gives; no cutoff)
  cutoff  a rank is called only when the best hit's ANI is >= GmG's per-rank ANI cutoff (0.7109 x3, 0.7715, 0.90, 0.9938)
GmG itself is read from `gmg run` results.tsv. Output: one row per (tool, variant, set, mode, group, rank): n, called, precision, recall.
ANI is each tool's own estimate (skani ANI, mash 1 - distance, sourmash query containment ANI, GmG containment ANI): comparable in size, not identical."""
import csv, glob, os, sys, collections

cdir, repo, out = sys.argv[1:4]
R = ["phylum", "class", "order", "family", "genus", "species"]
CUT = [0.7109, 0.7109, 0.7109, 0.7715, 0.90, 0.9938]
meta_p = sys.argv[4] if len(sys.argv) > 4 else None
if not meta_p:
    meta_p = glob.glob(os.path.dirname(cdir) + "/release/2026.10.2/meta_reps.tsv")[0]
meta = {r["id"]: r for r in csv.DictReader(open(meta_p), delimiter="\t")}
def tsv(p): return list(csv.DictReader(open(p), delimiter="\t"))
truth = {
    "s1": {r["id"]: r for r in tsv(f"{cdir}/sets/s1_post_freeze.tsv")},
    "s2": {r["id"]: r for r in tsv(f"{cdir}/sets/s2_truth.tsv")},
    "s3": {r["id"]: r for r in tsv(f"{cdir}/sets/s3_truth.tsv")},
}
stem = lambda p: os.path.basename(p).replace(".fa.gz", "").replace(".fna.gz", "")


def hits_skani(d):
    h = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{d}/hits.tsv"), delimiter="\t"):
        h[stem(r["Query_file"])].append((stem(r["Ref_file"]), float(r["ANI"]) / 100))
    return h


def hits_mash(d):
    h = collections.defaultdict(list)
    for line in open(f"{d}/hits.tsv"):
        a = line.rstrip("\n").split("\t")
        h[stem(a[1])].append((stem(a[0]), 1 - float(a[2])))
    return h


def hits_mashscreen(d):
    h = collections.defaultdict(list)
    for line in open(f"{d}/hits.tsv"):
        a = line.rstrip("\n").split("\t")  # query, identity, shared, multiplicity, p, ref-id
        h[a[0]].append((stem(a[5]), float(a[1])))
    return h


def hits_sourmash(d):
    h = collections.defaultdict(list)
    for f in glob.glob(f"{d}/gat/*.csv"):
        q = stem(f).replace(".csv", "")
        for i, r in enumerate(csv.DictReader(open(f))):
            ani = r.get("query_containment_ani") or r.get("average_containment_ani") or "0"
            name = r["name"].split()[0]
            h[os.path.basename(f)[:-4]].append((name, float(ani or 0)))
    return h


READERS = {"skani": hits_skani, "mash": hits_mash, "mashscreen": hits_mashscreen, "sourmash": hits_sourmash}


def exclude(q, qrow, mode, hit):
    """True if the query may not see this reference genome."""
    if mode == "none":
        return False
    if hit == qrow.get("bfd_id"):
        return True
    m = meta.get(hit)
    if m is None or mode == "loo":
        return False
    if mode == "los":
        return bool(qrow["species"]) and m["species"] == qrow["species"]
    if mode == "log":
        return bool(qrow["genus"]) and m["genus"] == qrow["genus"]
    return False


def calls_from(hits, variant):
    """best remaining hit -> (calls dict, ani)"""
    if not hits:
        return {r: "" for r in R}, 0.0
    hid, ani = max(hits, key=lambda x: x[1])
    m = meta.get(hid)
    if not m:
        return {r: "" for r in R}, ani
    return {r: (m[r] if variant == "raw" or ani >= CUT[i] else "") for i, r in enumerate(R)}, ani


rows = []  # (tool, variant, set, mode, group, rank, n, called, precision, recall)


def metrics(tool, variant, s, mode, group, items):
    for rk in R:
        t = [(c, tr) for c, tr in items if tr.get(rk)]
        if not t:
            continue
        called = [(c, tr) for c, tr in t if c.get(rk)]
        right = sum(c[rk] == tr[rk] for c, tr in called)
        rows.append((tool, variant, s, mode, group, rk, len(t), len(called) / len(t), right / len(called) if called else "", sum(c.get(rk) == tr[rk] for c, tr in t) / len(t)))


def group_of(s, tr):
    if s == "s1":
        return "species_in_bfd=" + tr["species_in_bfd"] if tr["species_in_bfd"] == "True" else "genus_in_bfd=" + tr["genus_in_bfd"]
    if s == "s3":
        return tr.get("variant", "") + (":" + tr.get("param", "") if tr.get("variant") in ("contam",) else "")
    return "all"


for tool in ["skani", "mash", "mashscreen", "sourmash"]:
    for s in ["s1", "s2", "s3"]:
        d = f"{cdir}/out/{tool}/{s}"
        if not os.path.exists(d):
            continue
        try:
            H = READERS[tool](d)
        except FileNotFoundError:
            continue
        modes = ["loo", "los", "log"] if s == "s2" else ["none"]
        for mode in modes:
            for variant in ("raw", "cutoff"):
                items = []
                for q, tr in truth[s].items():
                    if s == "s3" and tr.get("variant") == "negative":
                        continue
                    hits = [h for h in H.get(q, []) if not exclude(q, tr, mode, h[0])]
                    c, ani = calls_from(hits, variant)
                    items.append((c, tr))
                metrics(tool, variant, s, mode, "all", items)
                gs = collections.defaultdict(list)
                for c, tr in items:
                    gs[group_of(s, tr)].append((c, tr))
                if s != "s2":
                    for g, it in gs.items():
                        metrics(tool, variant, s, mode, g, it)
                if s == "s3":  # negative controls: any fungal call is a failure
                    neg = [q for q, tr in truth[s].items() if tr.get("variant") == "negative"]
                    n_called = 0
                    for q in neg:
                        c, ani = calls_from(H.get(q, []), variant)
                        n_called += bool(c["phylum"])
                    rows.append((tool, variant, s, mode, "negative_controls", "phylum", len(neg), n_called / len(neg) if neg else "", "", ""))

# GmG: results.tsv of `gmg run`
for s in ["s1", "s2", "s3"]:
    for mode, sub in ([("loo", "gmg.loo"), ("los", "gmg.los"), ("log", "gmg.log")] if s == "s2" else [("none", "gmg")]):
        p = f"{cdir}/out/gmg/{s}/{sub}/results.tsv"
        if not os.path.exists(p):
            continue
        res = {r["query_id"]: r for r in tsv(p)}
        items = []
        for q, tr in truth[s].items():
            if s == "s3" and tr.get("variant") == "negative":
                continue
            r = res.get(q)
            c = {rk: (r.get("call_" + rk, "") if r else "") for rk in R}
            items.append((c, tr))
        metrics("gmg", "native", s, mode, "all", items)
        if s != "s2":
            gs = collections.defaultdict(list)
            for c, tr in items:
                gs[group_of(s, tr)].append((c, tr))
            for g, it in gs.items():
                metrics("gmg", "native", s, mode, g, it)
        if s == "s3":
            neg = [q for q, tr in truth[s].items() if tr.get("variant") == "negative"]
            n_called = sum(bool(res.get(q, {}).get("call_phylum")) for q in neg)
            rows.append(("gmg", "native", s, mode, "negative_controls", "phylum", len(neg), n_called / len(neg) if neg else "", "", ""))

with open(out, "w") as o:
    w = csv.writer(o, delimiter="\t", lineterminator="\n")
    w.writerow(["tool", "variant", "set", "mode", "group", "rank", "n", "called", "precision", "recall"])
    for r in rows:
        w.writerow([f"{x:.4f}" if isinstance(x, float) else x for x in r])
print(len(rows), "rows ->", out)
