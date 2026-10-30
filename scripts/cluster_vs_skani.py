#!/usr/bin/env python3
"""Check the 95% ANI clusters of `gmg cluster` against skani (#47).

usage:
  cluster_vs_skani.py sample CLUSTERS.tsv GENOME_DIR OUT_DIR [--seed N] [--strata N] [--complex-cap N]
      writes OUT_DIR/lists/<name>.txt (one FASTA path per line) and OUT_DIR/manifest.tsv (kind, name, n, genus)
  cluster_vs_skani.py score CLUSTERS.tsv OUT_DIR REPORT.md [--prefix gmg_vs_skani]
      reads OUT_DIR/skani/<name>.tsv (`skani triangle -E` edge lists) and writes the report and OUT_DIR/summary.tsv

CLUSTERS.tsv = `gmg cluster` output: id, cluster (centroid id), is_centroid, ani_to_centroid, species, genus, n_hashes.
Three questions, all with skani (a different, alignment-based ANI) as the reference:
  1. members: is skani ANI(member, centroid) >= 95% where gmg's containment ANI says so? (stratified random clusters)
  2. chaining: inside one cluster, how many member pairs are < 95% (members each close to the centroid but not to each other)?
     (the largest clusters of species-complex genera)
  3. splits: do the centroids of different clusters of one genus reach >= 95% with each other (clusters that skani would merge)?
"""
import argparse, collections, csv, math, os, random, statistics, sys

IDENT = 0.95


def load(clusters):
    rows = list(csv.DictReader(open(clusters), delimiter="\t"))
    by = collections.defaultdict(list)
    for r in rows:
        by[r["cluster"]].append(r)
    return rows, by


def genome_path(gdir, gid):
    for suf in (".fa.gz", ".fna.gz", ".fa", ".fna", ".fasta.gz"):
        for sub in ("clean/", ""):
            p = f"{gdir}/{sub}{gid}{suf}"
            if os.path.exists(p):
                return p
    return None


def sample(a):
    rows, by = load(a.clusters)
    rng = random.Random(a.seed)
    os.makedirs(f"{a.out}/lists", exist_ok=True)
    man = []

    def write(kind, name, ids, genus=""):
        paths = [genome_path(a.genomes, i) for i in ids]
        ok = [(i, p) for i, p in zip(ids, paths) if p]
        if len(ok) < 2:
            return
        with open(f"{a.out}/lists/{name}.txt", "w") as f:
            f.write("\n".join(p for _, p in ok) + "\n")
        with open(f"{a.out}/lists/{name}.ids", "w") as f:
            f.write("\n".join(i for i, _ in ok) + "\n")
        man.append((kind, name, len(ok), genus))

    def members(cl, cap):
        ids = [r["id"] for r in by[cl]]
        rest = [i for i in ids if i != cl]
        rng.shuffle(rest)
        return [cl] + rest[: cap - 1]

    # 1. stratified random clusters of size >= 2
    strata = [(2, 2), (3, 5), (6, 20), (21, 100), (101, 10**9)]
    for lo, hi in strata:
        pool = [c for c, v in by.items() if lo <= len(v) <= hi]
        rng.shuffle(pool)
        for c in pool[: a.strata]:
            write("stratified", f"s{lo}_{c}", members(c, 30))
    # 2. the largest clusters of complex genera
    genera = ["Aspergillus", "Fusarium", "Saccharomyces", "Candida", "Candidozyma", "Cryptococcus", "Penicillium", "Trichophyton"]
    cen = {c: v for c, v in by.items()}
    for g in genera:
        cl = [c for c, v in cen.items() if sum(1 for x in v if x["genus"] == g) > len(v) / 2]
        cl.sort(key=lambda c: -len(cen[c]))
        for c in cl[:4]:
            if len(cen[c]) >= 5:
                write("complex", f"c_{g}_{c}", members(c, a.complex_cap), g)
        # 3. all centroids of the genus (capped): do separate clusters reach 95% with each other?
        cs = [c for c in cl]
        if len(cs) >= 2:
            rng.shuffle(cs)
            write("centroids", f"k_{g}", cs[:150], g)
    with open(f"{a.out}/manifest.tsv", "w") as f:
        f.write("kind\tname\tn\tgenus\n")
        for m in man:
            f.write("\t".join(map(str, m)) + "\n")
    print(len(man), "lists", collections.Counter(m[0] for m in man), file=sys.stderr)


def read_edges(path, ids_by_path):
    """skani -E output: Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name"""
    ani = {}
    if not os.path.exists(path):
        return None
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            a, b = ids_by_path.get(c[0]), ids_by_path.get(c[1])
            if a and b:
                v = float(c[2]) / 100.0
                ani[(a, b)] = ani[(b, a)] = v
    return ani


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def pct(x):
    return "-" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{100 * x:.1f}%"


def score(a):
    rows, by = load(a.clusters)
    meta = {r["id"]: r for r in rows}
    man = list(csv.DictReader(open(f"{a.out}/manifest.tsv"), delimiter="\t"))
    md = ["# GmG 95% clusters against skani (#47)\n"]
    summ = []
    pairs_all = []       # member vs centroid: (gmg ani, skani ani or None)
    per_cluster = []
    centroid_rows = []
    for m in man:
        ids = open(f"{a.out}/lists/{m['name']}.ids").read().split()
        paths = open(f"{a.out}/lists/{m['name']}.txt").read().split()
        e = read_edges(f"{a.out}/skani/{m['name']}.tsv", dict(zip(paths, ids)))
        if e is None:
            continue
        if m["kind"] == "centroids":
            near = collections.Counter()
            hits = []
            for (x, y), v in e.items():
                if x < y and v >= IDENT:
                    hits.append((x, y, v))
            n_pairs = len(ids) * (len(ids) - 1) // 2

            def ratio(x, y):
                a, b = int(meta[x]["n_hashes"]), int(meta[y]["n_hashes"])
                return max(a, b) / max(1, min(a, b))
            allr = sorted(ratio(x, y) for i, x in enumerate(ids) for y in ids[i + 1:])
            hr = sorted(ratio(x, y) for x, y, v in hits)
            centroid_rows.append((m["genus"], len(ids), n_pairs, hits, statistics.median(allr) if allr else float("nan"),
                                  statistics.median(hr) if hr else float("nan")))
            continue
        cl = meta[ids[0]]["cluster"]
        n = len(ids)
        below_c, tot_c = 0, 0
        for x in ids[1:]:
            g = float(meta[x]["ani_to_centroid"])
            s = e.get((ids[0], x))
            if m["kind"] == "stratified":
                pairs_all.append((g, s))
            tot_c += 1
            below_c += (s is None or s < IDENT)
        pw = [(e.get((x, y)) or 0.0) for i, x in enumerate(ids) for y in ids[i + 1:]]
        below = sum(1 for v in pw if v < IDENT)
        nh = sorted(int(meta[i]["n_hashes"]) for i in ids[1:])
        size_ratio = int(meta[ids[0]]["n_hashes"]) / nh[len(nh) // 2] if nh and nh[len(nh) // 2] else float("nan")
        per_cluster.append(dict(size_ratio=size_ratio, kind=m["kind"], name=m["name"], genus=m["genus"], cluster=cl, n=n, centroid_below=below_c,
                                centroid_total=tot_c, pairs=len(pw), pairs_below=below, min_pair=min(pw) if pw else None,
                                species=len({meta[i]["species"] for i in ids}), species_labels="; ".join(sorted({meta[i]["species"] for i in ids})[:6])))
    # 1. members vs centroid
    strat = [c for c in per_cluster if c["kind"] == "stratified"]
    k = sum(c["centroid_below"] for c in strat)
    nt = sum(c["centroid_total"] for c in strat)
    lo, hi = wilson(k, nt)
    md.append(f"## 1. Members against their centroid (stratified random clusters)\n\n{len(strat)} clusters (sizes 2 to >100, up to 30 members each), "
              f"{nt} member-centroid pairs. gmg's containment ANI (k=21, scaled 10000) is >= 0.95 for all of them by construction. **skani ANI "
              f"< 95% (or not reported, i.e. < ~80% / too little alignment) for {k} of {nt} ({pct(k / nt if nt else None)}; 95% CI {pct(lo)} to {pct(hi)}).**\n")
    both = [(g, s) for g, s in pairs_all if s is not None]
    if both:
        d = [s - g for g, s in both]
        md.append(f"Where skani reports a value ({len(both)} pairs): mean skani - gmg = {statistics.mean(d):+.4f}, median {statistics.median(d):+.4f}, "
                  f"5th percentile {sorted(d)[int(0.05 * len(d))]:+.4f}, 95th {sorted(d)[int(0.95 * len(d))]:+.4f}. "
                  f"skani ANI < 0.95 for {sum(1 for g, s in both if s < IDENT)} of them.\n")
        bins = [(0.95, 0.96), (0.96, 0.97), (0.97, 0.98), (0.98, 0.99), (0.99, 1.01)]
        md.append("| gmg ANI to centroid | n | mean skani ANI | skani < 0.95 |\n|---|---|---|---|")
        for lo_, hi_ in bins:
            s = [(g, sk) for g, sk in both if lo_ <= g < hi_]
            if s:
                md.append(f"| [{lo_:.2f}, {hi_:.2f}) | {len(s)} | {statistics.mean(x[1] for x in s):.4f} | {sum(1 for x in s if x[1] < IDENT)} |")
        md.append("")
    # 2. chaining in species-complex clusters
    cx = [c for c in per_cluster if c["kind"] == "complex"]
    if cx:
        md.append("## 2. Chaining inside the largest clusters of species-complex genera\n\nAll pairs of up to 60 members (centroid included). "
                  "`pairs < 95%` counts member pairs that skani puts below 95% although both are within 95% (gmg) of the centroid.\n")
        md.append("| genus | cluster (centroid) | size in index | members tested | member-centroid < 95% | pairs < 95% | min pair ANI (0 = not reported) | centroid size / median member | species labels |\n|---|---|---|---|---|---|---|---|---|")
        for c in sorted(cx, key=lambda c: (c["genus"], -c["n"])):
            minp = ("-" if c["min_pair"] is None else "%.3f" % c["min_pair"]) + f" | {c['size_ratio']:.2f}"
            md.append(f"| {c['genus']} | {c['cluster'][:34]} | {len(by[c['cluster']])} | {c['n']} | {c['centroid_below']}/{c['centroid_total']} | "
                      f"{c['pairs_below']}/{c['pairs']} ({pct(c['pairs_below'] / c['pairs'] if c['pairs'] else None)}) | "
                      f"{minp} | {c['species']}: {c['species_labels'][:90]} |")
        md.append("")
    # 3. splits
    if centroid_rows:
        md.append("## 3. Do separate clusters of one genus merge at 95% under skani?\n\nAll centroid pairs of the clusters of a genus (up to 150 centroids). "
                  "gmg's greedy clustering keeps centroids below 95% containment ANI of each other; skani may disagree.\n")
        md.append("| genus | centroids | pairs | skani ANI >= 95% | >= 97% | size ratio of those pairs (median, larger/smaller n_hashes) | same for all pairs | example (ANI) |\n|---|---|---|---|---|---|---|---|")
        for g, n, npairs, hits, r_all, r_hit in centroid_rows:
            ex = "; ".join(f"{meta[x]['species']} / {meta[y]['species']} ({v:.3f})" for x, y, v in sorted(hits, key=lambda t: -t[2])[:2])
            md.append(f"| {g} | {n} | {npairs} | {len(hits)} | {sum(1 for h in hits if h[2] >= 0.97)} | {'-' if r_hit != r_hit else '%.2f' % r_hit} | {r_all:.2f} | {ex[:110]} |")
        md.append("")
    with open(f"{a.out}/summary.tsv", "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["kind", "name", "genus", "n", "centroid_below", "centroid_total", "pairs", "pairs_below", "min_pair", "species_labels"])
        for c in per_cluster:
            w.writerow([c["kind"], c["name"], c["genus"], c["n"], c["centroid_below"], c["centroid_total"], c["pairs"], c["pairs_below"],
                        "" if c["min_pair"] is None else f"{c['min_pair']:.4f}", c["species_labels"]])
    open(a.report, "w").write("\n".join(md) + "\n")
    print("\n".join(md))


ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)
s = sub.add_parser("sample"); s.add_argument("clusters"); s.add_argument("genomes"); s.add_argument("out")
s.add_argument("--seed", type=int, default=5); s.add_argument("--strata", type=int, default=25); s.add_argument("--complex-cap", type=int, default=60)
t = sub.add_parser("score"); t.add_argument("clusters"); t.add_argument("out"); t.add_argument("report")
a = ap.parse_args()
sample(a) if a.cmd == "sample" else score(a)
