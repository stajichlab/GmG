#!/usr/bin/env python3
"""Score the chimera/hybrid, small-input and negative-control runs of the second stress test (#41, #45).

usage: stress2_score.py chim RESULTS.tsv HITS.tsv TRUTH.tsv OUT.md
       stress2_score.py small RESULTS.tsv SMALL.tsv OUT.md
       stress2_score.py neg RESULTS.tsv NEGATIVES.txt OUT.md
Per-genome output stays on the cluster; OUT.md holds aggregates only (plus the small-input and negative rows, which are one line per organism).
"""
import collections, csv, statistics, sys

L = ["phylum", "class", "order", "family", "genus", "species"]


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def chim(res_p, hits_p, truth_p, out):
    res = {r["query_id"]: r for r in tsv(res_p)}
    hits = collections.defaultdict(list)
    for h in tsv(hits_p):
        hits[h["query_id"]].append(h)
    md = ["## Chimeras and hybrids\n"]
    groups = collections.defaultdict(list)
    for t in tsv(truth_p):
        if t["id"] in res:
            groups[(t["kind"], t["variant"], int(t["pct"]))].append(t)
    md.append("Per group: n; status; the call at the last rank the two parents share (correct / called); at the rank where they differ the call is the "
              "majority parent's taxon, the minority's, another taxon, or none; whether the minority parent's taxon appears among the top hits "
              "(`hits.tsv`, top 10) and how often the majority's does.\n")
    md.append("| kind | variant | % A | n | status ok / low / no_match | common rank | called (correct) | at split: majority / minority / other / none | minority in top hits | majority in top hits | median ANI | flags |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for key in sorted(groups, key=lambda k: (k[0], k[1] == "hybrid", -k[2])):
        kind, var, pct = key
        g = groups[key]
        sr = {"same_genus": "species", "same_order": "genus", "same_class": "order"}[kind]
        common = L[max(0, L.index(sr) - 1)]
        st = collections.Counter(res[t["id"]]["status"] for t in g)
        called = sum(bool(res[t["id"]]["call_" + common]) for t in g)
        right = sum(res[t["id"]]["call_" + common] == t[common] and bool(res[t["id"]]["call_" + common]) for t in g)
        maj = mino = oth = none = 0
        inmin = inmaj = 0
        for t in g:
            r = res[t["id"]]
            c = r["call_" + sr]
            a, b = t[sr], t["b_" + sr]
            if not c: none += 1
            elif c == a and c != b: maj += 1
            elif c == b and c != a: mino += 1
            elif c == a: maj += 1
            else: oth += 1
            hs = hits[t["id"]]
            if b and any(h[sr] == b for h in hs) and b != a: inmin += 1
            if a and any(h[sr] == a for h in hs): inmaj += 1
        anis = sorted(float(res[t["id"]]["ani"]) for t in g if res[t["id"]]["ani"])
        fl = collections.Counter(f for t in g for f in res[t["id"]]["qc_flags"].split(",") if f)
        md.append(f"| {kind} | {var} | {pct} | {len(g)} | {st['ok']} / {st['low_confidence']} / {st['no_match']} | {common} | {called} ({right}) | "
                  f"{maj} / {mino} / {oth} / {none} (rank: {sr}) | {inmin}/{len(g) if sr != 'species' else 'n/a'} | {inmaj}/{len(g)} | "
                  f"{statistics.median(anis) if anis else float('nan'):.4f} | {', '.join(f'{k} {v}' for k, v in fl.most_common(3))} |")
    open(out, "w").write("\n".join(md) + "\n")
    print("\n".join(md))


def small(res_p, small_p, out):
    res = {r["query_id"]: r for r in tsv(res_p)}
    md = ["## Mitochondrion-, plasmid- and virus-only inputs\n",
          "| input | kind | length | status | call | n_hashes | flags |", "|---|---|---|---|---|---|---|"]
    for s in tsv(small_p):
        r = res.get(s["id"])
        if not r:
            md.append(f"| {s['title'][:60]} | {s['kind']} | {s['length']} | not run |  |  |  |"); continue
        call = next((r["call_" + x] for x in reversed(L) if r["call_" + x]), "") or ""
        md.append(f"| {s['title'][:60]} | {s['kind']} | {int(s['length']):,} | {r['status']} | {call} | {r['n_hashes']} | {r['qc_flags']} |")
    open(out, "w").write("\n".join(md) + "\n")
    print("\n".join(md))


def neg(res_p, neg_p, out):
    res = {r["query_id"]: r for r in tsv(res_p)}
    md = ["## Archaea and more non-fungal negatives\n",
          "| genome | status | deepest call | ANI | flags | verdict |", "|---|---|---|---|---|---|"]
    bad = 0
    for line in open(neg_p):
        if line.startswith("#") or not line.strip(): continue
        acc, name = line.split(None, 1)
        r = res.get(acc)
        if not r:
            md.append(f"| {name.strip()} | not run |  |  |  |  |"); continue
        call = next((r["call_" + x] for x in reversed(L) if r["call_" + x]), "") or ""
        flagged = any(f in r["qc_flags"] for f in ("not_fungal", "distant_non_fungal_hit"))
        confident = bool(r["call_genus"] or r["call_species"]) and not flagged
        verdict = "no_match" if r["status"] == "no_match" else "flagged" if flagged else "WEAK fungal phylum call" if call and not confident else "CONFIDENT fungal call"
        bad += verdict == "CONFIDENT fungal call"
        md.append(f"| {name.strip()[:60]} | {r['status']} | {call} | {r['ani'][:6]} | {r['qc_flags']} | {verdict} |")
    md.append(f"\nConfident fungal calls: {bad}.")
    open(out, "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    mode = sys.argv[1]
    {"chim": lambda: chim(*sys.argv[2:6]), "small": lambda: small(*sys.argv[2:5]), "neg": lambda: neg(*sys.argv[2:5])}[mode]()
