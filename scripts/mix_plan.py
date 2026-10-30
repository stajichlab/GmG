#!/usr/bin/env python3
"""Plan the mixed-read experiment (#45): which real read sets are mixed with which, at which distance.

usage: mix_plan.py TITR_DIR COCCI_TRUTH COCCI_DIR OUT_DIR [LOW_PAIRS=600000] [HIGH_PAIRS=2700000] [SEED=11]
  TITR_DIR   the reads-titration run dir (truth.tsv, labels.tsv, sheet.tsv, n1500000/results.tsv): ZyGoLife read sets, truth = the lab's
             file-name genus resolved through the taxdump; only sets whose genus is called correctly alone at 1.5M pairs are used
  COCCI_TRUTH data/validation/cocci_reads.truth.tsv (Coccidioides read sets with a curated species); COCCI_DIR the folder of their FASTQ
Writes OUT_DIR/plan.tsv (one row per pair: distance same_species | same_genus | same_order | diff_phylum, parents A and B with their
lineage and read files, `high` = both parents have >= HIGH_PAIRS read pairs so the high-depth variants can be made) and OUT_DIR/triples.tsv
(three parents each). Read-pair counts are cached in OUT_DIR/pairs.tsv (zcat | wc -l, slow: run it in a batch job). Read-only on the inputs.
"""
import collections, csv, os, random, subprocess, sys

titr, cocci_truth, cocci_dir, out = sys.argv[1:5]
LOW = int(sys.argv[5]) if len(sys.argv) > 5 else 600_000
HIGH = int(sys.argv[6]) if len(sys.argv) > 6 else 2_700_000
rng = random.Random(int(sys.argv[7]) if len(sys.argv) > 7 else 11)
os.makedirs(out, exist_ok=True)
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))

cache_p = f"{out}/pairs.tsv"
cache = {r["path"]: int(r["pairs"]) for r in rd(cache_p)} if os.path.exists(cache_p) else {}


def count_pairs(path):
    if path not in cache:
        n = subprocess.run(f"zcat '{path}' | wc -l", shell=True, capture_output=True, text=True).stdout.strip()
        cache[path] = int(n) // 4
        with open(cache_p, "w") as f:
            f.write("path\tpairs\n" + "".join(f"{k}\t{v}\n" for k, v in cache.items()))
    return cache[path]


hosts = []  # dicts: id, r1, r2, lineage, label (species name), curated (species is curated)
truth = {r["id"]: r for r in rd(f"{titr}/truth.tsv")}
label = {r["sample"]: r["name"] for r in rd(f"{titr}/labels.tsv")}
files = collections.defaultdict(list)
for r in rd(f"{titr}/sheet.tsv"):
    files[r["id"]].append(r["path"])
res = {r["query_id"]: r for r in rd(f"{titr}/n1500000/results.tsv")}
for i, t in truth.items():
    r = res.get(i)
    if not r or r["call_genus"] != t["genus"] or len(files[i]) < 2:
        continue
    r1, r2 = sorted(files[i])[:2]
    name = label[i]
    hosts.append({"id": i, "r1": r1, "r2": r2, "lin": {k: t[k] for k in RANKS[:5]} | {"species": name if not name.endswith(" sp.") else ""},
                  "curated": False, "name": name})
cocci = [t for t in rd(cocci_truth) if int(t["bytes"]) >= 900_000_000]
rng.shuffle(cocci)
per_species = collections.Counter()
for t in cocci:
    if per_species[t["species"]] >= 10:     # counting read pairs decompresses the whole file: a few candidates per species are enough
        continue
    r1 = f"{cocci_dir}/{t['id']}.fastq.gz"
    r2 = r1.replace("_R1_", "_R2_")
    if os.path.exists(r1) and os.path.exists(r2):
        per_species[t["species"]] += 1
        hosts.append({"id": t["id"], "r1": r1, "r2": r2, "lin": {k: t[k] for k in RANKS}, "curated": True, "name": t["species"]})
for h in hosts:
    h["pairs"] = count_pairs(h["r1"])       # R1 only (R2 has the same number of pairs; mix_reads.py fails loudly if it does not)
hosts = [h for h in hosts if h["pairs"] >= LOW + LOW // 6]
rng.shuffle(hosts)
print(len(hosts), "hosts with the correct genus alone and >=", LOW, "read pairs", file=sys.stderr)

used = collections.Counter()
chosen = []   # (distance, A, B)


def pick(distance, cands, limit):
    n = 0
    for a, b in cands:
        if n >= limit or used[a["id"]] >= 2 or used[b["id"]] >= 2:
            continue
        used[a["id"]] += 1; used[b["id"]] += 1
        chosen.append((distance, a, b)); n += 1


named = collections.defaultdict(list)       # species name -> hosts (named, not "sp.")
for h in hosts:
    if h["lin"]["species"]:
        named[h["lin"]["species"]].append(h)
same_sp = [(v[0], v[1]) for k, v in sorted(named.items()) if len(v) >= 2]
pick("same_species", same_sp, 8)
genus_sp = collections.defaultdict(dict)
for h in hosts:
    if h["lin"]["species"]:
        genus_sp[h["lin"]["genus"]].setdefault(h["lin"]["species"], h)
same_ge = []
for g, d in sorted(genus_sp.items()):
    sp = sorted(d)
    if len(sp) >= 2:
        same_ge.append((d[sp[0]], d[sp[1]]))
pick("same_genus", same_ge, 8)
pool = hosts[:]
same_or = [(a, b) for a in pool for b in pool if a["id"] < b["id"] and a["lin"]["order"] == b["lin"]["order"] and a["lin"]["genus"] != b["lin"]["genus"]]
rng.shuffle(same_or)
pick("same_order", same_or, 8)
diff_ph = [(a, b) for a in pool for b in pool if a["id"] < b["id"] and a["lin"]["phylum"] != b["lin"]["phylum"]]
rng.shuffle(diff_ph)
pick("diff_phylum", diff_ph, 8)

cols = ["pair_id", "distance", "high", "a_id", "a_r1", "a_r2", "b_id", "b_r1", "b_r2"] + [f"a_{r}" for r in RANKS] + [f"b_{r}" for r in RANKS]
with open(f"{out}/plan.tsv", "w") as f:
    f.write("\t".join(cols) + "\n")
    for n, (d, a, b) in enumerate(chosen):
        high = int(a["pairs"] >= HIGH + HIGH // 6 and b["pairs"] >= HIGH + HIGH // 6)
        row = [f"p{n:02d}", d, high, a["id"], a["r1"], a["r2"], b["id"], b["r1"], b["r2"]] + [a["lin"][r] for r in RANKS] + [b["lin"][r] for r in RANKS]
        f.write("\t".join(map(str, row)) + "\n")
print(collections.Counter(d for d, a, b in chosen), "pairs;", sum(1 for d, a, b in chosen if a["pairs"] >= HIGH and b["pairs"] >= HIGH), "with high depth", file=sys.stderr)

# triples: three genera of one order; two species of a genus + another order of the phylum; three phyla
trip = []   # (type, (host, host, host))
usedt = set()


def free(h):
    return h["id"] not in usedt


def take(ty, *hs):
    usedt.update(h["id"] for h in hs)
    trip.append((ty, hs))


for _ in range(2):   # three genera of one order
    byo = collections.defaultdict(dict)
    for h in hosts:
        if free(h):
            byo[(h["lin"]["phylum"], h["lin"]["order"])].setdefault(h["lin"]["genus"], h)
    c = [list(v.values())[:3] for v in byo.values() if len(v) >= 3]
    if c:
        take("same_order", *c[0])
for _ in range(2):   # two species of one genus plus a host of another order of the phylum
    g2 = [(a, b) for a, b in same_ge if free(a) and free(b)]
    if g2:
        a, b = g2[0]
        o = [h for h in hosts if free(h) and h["lin"]["phylum"] == a["lin"]["phylum"] and h["lin"]["order"] != a["lin"]["order"]]
        if o:
            take("genus_pair_plus_order", a, b, o[0])
for _ in range(2):   # three phyla
    byp = {}
    for h in hosts:
        if free(h):
            byp.setdefault(h["lin"]["phylum"], h)
    if len(byp) >= 3:
        take("three_phyla", *list(byp.values())[:3])
with open(f"{out}/triples.tsv", "w") as f:
    f.write("\t".join(["triple_id", "type"] + [f"{p}_{k}" for p in "abc" for k in ["id", "r1", "r2"] + RANKS]) + "\n")
    for n, (ty, hs) in enumerate(trip):
        row = [f"t{n:02d}", ty]
        for h in hs:
            row += [h["id"], h["r1"], h["r2"]] + [h["lin"][r] for r in RANKS]
        f.write("\t".join(map(str, row)) + "\n")
print(len(trip), "triples", file=sys.stderr)
