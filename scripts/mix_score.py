#!/usr/bin/env python3
"""Score the mixed-read experiment (#45): what the call, the ANI and the hit list say about a two- or three-fungus read set.

usage: mix_score.py PLAN.tsv TRIPLES.tsv --run LABEL=RESULTS.tsv:HITS.tsv [--run ...]    (markdown on stdout)
Variant ids are PAIR__DEPTH__TAG: TAG is A<pct> (pct of the read pairs from parent A: A90 = 90:10, A100 = A alone, A0 = B alone) for pairs and
`mix` for triples (shares in the triples table). Truth = the parents' lineages in the plan (the lab's file-name genus resolved through the
taxdump; species only for the Coccidioides read sets, which have a curated species).
Per distance (same_species / same_genus / same_order / diff_phylum), depth and share of parent A, at the DECISIVE rank (the first rank where the two
parents differ; genus when they share it): the call is the majority parent's taxon, the minority parent's, another taxon, or absent; the minority taxon
is in the top-10 hits or not; the top-hit ANI relative to A alone against p^(1/21) (p = A's share); retention of A's genus; the `weak_phylum` flag.
"""
import collections, csv, statistics, sys

RANKS = ["phylum", "class", "order", "family", "genus", "species"]
args = sys.argv[1:]
plan = {r["pair_id"]: r for r in csv.DictReader(open(args[0]), delimiter="\t")}
trip = {r["triple_id"]: r for r in csv.DictReader(open(args[1]), delimiter="\t")}
runs = []
for i, x in enumerate(args):
    if x == "--run":
        label, rest = args[i + 1].split("=", 1)
        rp, hp = rest.split(":")
        res = {r["query_id"]: r for r in csv.DictReader(open(rp), delimiter="\t")}
        hits = collections.defaultdict(list)
        for h in csv.DictReader(open(hp), delimiter="\t"):
            hits[h["query_id"]].append(h)
        for v in hits.values():
            v.sort(key=lambda h: int(h["hit_rank"]))
        runs.append((label, res, hits))


def lin(row, who):
    return {r: row[f"{who}_{r}"] for r in RANKS}


def decisive(pr):
    A, B = lin(pr, "a"), lin(pr, "b")
    for r in RANKS[:5]:
        if A[r] != B[r]:
            return r
    return "genus"


def outcome(call, A, B, r, pa):
    """majority | minority | shared | other | none for one rank."""
    if not call:
        return "none"
    if A[r] == B[r]:
        return "shared" if call == A[r] else "other"
    if pa == 50:   # a tie: either parent is a correct answer
        return "majority" if call in (A[r], B[r]) else "other"
    maj, mino = (A[r], B[r]) if pa > 50 else (B[r], A[r])
    return "majority" if call == maj else "minority" if call == mino else "other"


def pct(n, d):
    return "-" if not d else f"{100 * n / d:.0f}%"


out = []
for label, res, hits in runs:
    rows = []   # (plan row, depth, share of A, result, hits)
    for vid, r in res.items():
        parts = vid.split("__")
        if len(parts) == 3 and parts[0] in plan and parts[2].startswith("A"):
            rows.append((plan[parts[0]], parts[1], float(parts[2][1:]), r, hits.get(vid, [])))
    out.append(f"\n## Run `{label}`: {len(rows)} two-fungus inputs\n")
    pure = {(pr["pair_id"], d): r for pr, d, p, r, h in rows if p == 100.0}
    for depth in sorted({d for _, d, _, _, _ in rows}):
        out.append(f"\n### Depth `{depth}`\n")
        out.append("| distance | rank | A share | n | majority / shared | minority | other | no call | minority in top-10 hits | ANI vs A alone: observed / p^(1/21) | genus of A kept | weak_phylum |")
        out.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for dist in ["same_species", "same_genus", "same_order", "diff_phylum"]:
            shares = sorted({p for pr, d, p, r, h in rows if d == depth and pr["distance"] == dist}, reverse=True)
            for p in shares:
                L = [(pr, r, h) for pr, d, pp, r, h in rows if d == depth and pr["distance"] == dist and pp == p]
                dr = decisive(L[0][0])
                oc = collections.Counter()
                mi = mn = gk = gn = wk = 0
                ratios = []
                for pr, r, h in L:
                    A, B = lin(pr, "a"), lin(pr, "b")
                    oc[outcome(r.get(f"call_{dr}", ""), A, B, dr, p)] += 1
                    if A[dr] != B[dr] and 0 < p < 100:
                        mn += 1
                        mino = B[dr] if p >= 50 else A[dr]
                        mi += any(x.get(dr, "") == mino for x in h[:10])
                    base = pure.get((pr["pair_id"], depth))
                    if base is not None and base.get("ani") and r.get("ani") and p >= 50:
                        ratios.append((float(r["ani"]) / float(base["ani"]), (p / 100) ** (1 / 21)))
                    if base is not None and base.get("call_genus") == A["genus"]:
                        gn += 1
                        gk += r.get("call_genus", "") == A["genus"]
                    wk += "weak_phylum" in r.get("qc_flags", "")
                n = len(L)
                ani = "-" if not ratios else f"{statistics.mean(o for o, e in ratios):.4f} / {statistics.mean(e for o, e in ratios):.4f}"
                out.append(f"| {dist} | {dr} | {p:g}% | {n} | {pct(oc['majority'] + oc['shared'], n)} | {pct(oc['minority'], n)} | {pct(oc['other'], n)} | {pct(oc['none'], n)} | {pct(mi, mn)} | {ani} | {pct(gk, gn)} | {wk} |")
    cocci = [(pr, d, p, r) for pr, d, p, r, h in rows if pr["a_species"].startswith("Coccidioides") and pr["b_species"].startswith("Coccidioides")]
    if cocci:
        out.append("\n**Coccidioides pairs (curated species): the species called**\n")
        out.append("| distance | depth | A share | n | species of A | species of B | other or none |")
        out.append("|---|---|---|---|---|---|---|")
        for dist in ("same_species", "same_genus"):
            for depth in sorted({d for pr, d, p, r in cocci}):
                for p in sorted({p for pr, d, p, r in cocci if pr["distance"] == dist and d == depth}, reverse=True):
                    L = [(pr, r) for pr, d, pp, r in cocci if pr["distance"] == dist and d == depth and pp == p]
                    na = sum(r.get("call_species", "") == pr["a_species"] for pr, r in L)
                    nb = sum(r.get("call_species", "") == pr["b_species"] and pr["a_species"] != pr["b_species"] for pr, r in L)
                    out.append(f"| {dist} | {depth} | {p:g}% | {len(L)} | {pct(na, len(L))} | {pct(nb, len(L))} | {pct(len(L) - na - nb, len(L))} |")
    tr = [(trip[vid.split("__")[0]], r, hits.get(vid, [])) for vid, r in res.items() if vid.split("__")[0] in trip and vid.endswith("__mix")]
    if tr:
        out.append("\n**Three-genome mixtures (60:30:10)**\n")
        out.append("| triple | type | A (60%) | B (30%) | C (10%) | deepest call | genus call is A / B / C | genus of A / B / C in the top-10 hits |")
        out.append("|---|---|---|---|---|---|---|---|")
        for t, r, h in sorted(tr, key=lambda x: x[0]["triple_id"]):
            g = [t[f"{p}_genus"] for p in "abc"]
            deepest = next((r[f"call_{k}"] for k in reversed(RANKS) if r.get(f"call_{k}")), "-")
            hg = ["yes" if any(x.get("genus", "") == gg for x in h[:10]) else "no" for gg in g]
            gc = ["yes" if r.get("call_genus", "") == gg else "no" for gg in g]
            out.append(f"| {t['triple_id']} | {t['type']} | {g[0]} | {g[1]} | {g[2]} | {deepest} | {'/'.join(gc)} | {'/'.join(hg)} |")
print("\n".join(out))
