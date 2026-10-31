#!/usr/bin/env python3
"""Evaluate the mixture detector (#31) from results.jsonl files that carry qc.mixture (schema 8), at any threshold.

usage: mixture_eval.py ROOT [--f F] [--h H] [--cov C] [--rest R] [--sweep]
  ROOT = GmG_runs/mixdet/out with: asm/ (results.jsonl, truth.tsv), reads/<pair>/results.jsonl, truth_loo|los|log/*/results.jsonl, contam/out/results.jsonl
  --f/--h/--cov/--rest: the rule  second.explained_fraction >= F and second.hashes >= H and second.concentration >= C and second.fraction_of_rest >= R
Prints markdown tables: second-taxon fraction on the negatives, false-mixture rate per negative set, detection rate vs share per pair class,
abundance error of the explained fractions, and (--sweep) the false-mixture rate / detection of a grid of rules.
"""
import argparse, csv, glob, json, os, statistics, sys

ap = argparse.ArgumentParser()
ap.add_argument("root")
ap.add_argument("--f", type=float, default=0.01)
ap.add_argument("--h", type=int, default=20)
ap.add_argument("--cov", type=float, default=0.0)
ap.add_argument("--rest", type=float, default=0.0)
ap.add_argument("--sweep", action="store_true")
a = ap.parse_args()
R = a.root


def jl(pattern):
    out = []
    for p in sorted(glob.glob(pattern)):
        for line in open(p):
            if line.strip():
                out.append(json.loads(line))
    return out


def feat(r):
    """(explained, second fraction, hashes, ref_coverage, fraction_of_rest, second phylum, primary phylum-less) or None"""
    m = (r.get("qc") or {}).get("mixture")
    if not m:
        return None
    s = m.get("second")
    e = m["explained_fraction"]
    if not s:
        return dict(e=e, f=0.0, h=0, cov=0.0, conc=0.0, rest=0.0, taxon="", phylum="", n=r["qc"]["n_hashes"])
    return dict(e=e, f=s["explained_fraction"], h=s["hashes"], cov=s["ref_coverage"], conc=s.get("concentration", 0.0), rest=s["fraction_of_rest"],
                taxon=s["taxon"], phylum=s["phylum"], n=r["qc"]["n_hashes"])


def rule(x, f=None, h=None, cov=None, rest=None):
    f = a.f if f is None else f; h = a.h if h is None else h
    cov = a.cov if cov is None else cov; rest = a.rest if rest is None else rest
    return x is not None and x["f"] >= f and x["h"] >= h and x["conc"] >= cov and x["rest"] >= rest and x["f"] > 0


def pct(k, n):
    return f"{100 * k / n:.1f}%" if n else "n/a"


# ---------- sets
neg = {}   # name -> list of (id, feat)
for mode in ("loo", "los", "log"):
    rs = jl(f"{R}/truth_{mode}/*/results.jsonl")
    neg[f"truth set, {mode} (assemblies)"] = [(r["id"], feat(r)) for r in rs if r["status"] != "error"]
# further pure sets: ROOT/neg_<name>/**/results.jsonl (e.g. post-freeze assemblies that are not in the reference)
for d in sorted(glob.glob(f"{R}/neg_*")):
    rs = jl(f"{d}/results.jsonl") + jl(f"{d}/*/results.jsonl")
    neg[f"{os.path.basename(d)[4:]} (assemblies)"] = [(r["id"], feat(r)) for r in rs if r["status"] != "error"]
ct = {}
tp = f"{R}/contam/truth.tsv"
if os.path.exists(tp):
    for t in csv.DictReader(open(tp), delimiter="\t"):
        ct[t["id"]] = t
rs = jl(f"{R}/contam/out/results.jsonl")
for r in rs:
    t = ct.get(r["id"], {})
    k = "contaminated with " + t.get("donor", "?") if t.get("frac", "0") not in ("0", "0.0", "") else "hosts alone"
    neg.setdefault(f"contam sim: {k} (assemblies)", []).append((r["id"], feat(r)))

asm = {}
tp = f"{R}/asm/truth.tsv"
if os.path.exists(tp):
    for t in csv.DictReader(open(tp), delimiter="\t"):
        asm[t["id"]] = t
asm_rows = []
for r in jl(f"{R}/asm/results.jsonl"):
    t = asm.get(r["id"])
    if t and r["status"] != "error":
        asm_rows.append((t, feat(r), r))
for t, x, r in asm_rows:
    if t["b_pct"] in ("0", "100"):
        neg.setdefault("parents alone (assemblies, parents hidden)", []).append((t["id"], x))

plan = {}
pp = "/bigdata/stajichlab/shared/projects/GmG_runs/reads/mix/plan"
for fn in (f"{R}/plan.tsv", f"{R}/../../reads/mix/plan/plan.tsv", f"{pp}/plan.tsv"):
    if os.path.exists(fn):
        for t in csv.DictReader(open(fn), delimiter="\t"):
            plan[t["pair_id"]] = t
        break
reads_rows = []
for r in jl(f"{R}/reads/*/results.jsonl"):
    if r["status"] == "error":
        continue
    parts = r["id"].split("__")
    pid = parts[0]
    if len(parts) < 3 or pid not in plan:
        continue
    depth = parts[1]
    share = parts[2][1:]
    if share == "ix":      # triples
        continue
    reads_rows.append((plan[pid]["distance"], depth, float(share), feat(r), r))
for dist, depth, pa, x, r in reads_rows:
    if pa in (100.0, 0.0):
        neg.setdefault(f"pure read sets ({depth} depth)", []).append((r["id"], x))

# ---------- tables
print("## Negatives: second-taxon fraction (of all hashes) and hashes, per set\n")
print("| set | n | second fraction: median / 95th / 99th / max | second hashes: 99th / max | concentration of the second: 99th |")
print("|---|---|---|---|---|")
for name, v in neg.items():
    xs = [x for _, x in v if x]
    if not xs:
        continue
    q = lambda key, p: sorted(y[key] for y in xs)[min(len(xs) - 1, int(p * len(xs)))]
    print(f"| {name} | {len(xs)} | {q('f',.5):.4f} / {q('f',.95):.4f} / {q('f',.99):.4f} / {max(y['f'] for y in xs):.4f} | "
          f"{q('h',.99)} / {max(y['h'] for y in xs)} | {q('conc',.99):.3f} |")

print(f"\n## False-mixture rate at the rule (fraction >= {a.f}, hashes >= {a.h}, concentration >= {a.cov}, of rest >= {a.rest})\n")
print("| negative set | inputs | flagged | rate |")
print("|---|---|---|---|")
tot = [0, 0]
for name, v in neg.items():
    xs = [x for _, x in v if x]
    k = sum(rule(x) for x in xs)
    print(f"| {name} | {len(xs)} | {k} | {pct(k, len(xs))} |")

print("\n## Detection rate by minority share (assemblies, parents hidden); fraction of inputs flagged\n")
shares = sorted({float(t["b_pct"]) for t, _, _ in asm_rows if t["b_pct"] not in ("0", "100")})
print("| pair class | n pairs | " + " | ".join(f"{s:g}%" for s in shares) + " |")
print("|---|---|" + "---|" * len(shares))
for kind in ("same_genus", "same_order", "diff_phylum"):
    cells = []
    npairs = len({t["pair"] for t, _, _ in asm_rows if t["kind"] == kind})
    for s in shares:
        v = [x for t, x, _ in asm_rows if t["kind"] == kind and float(t["b_pct"]) == s]
        cells.append(f"{sum(rule(x) for x in v)}/{len(v)}")
    print(f"| {kind} | {npairs} | " + " | ".join(cells) + " |")

# a minority the reference does not cover cannot be explained: stratify by how well each parent alone is explained
alone = {}
for t, x, _ in asm_rows:
    if x and t["b_pct"] in ("0", "100"):
        alone[(t["pair"], t["b_pct"])] = x["e"]
KNOWN = 0.5
print(f"\n### the same, only pairs whose two parents are each explained >= {KNOWN:g} by their genus alone (the genus is well represented in the reference)\n")
print("| pair class | n pairs (of all) | " + " | ".join(f"{s:g}%" for s in shares) + " |")
print("|---|---|" + "---|" * len(shares))
for kind in ("same_genus", "same_order", "diff_phylum"):
    allp = {t["pair"] for t, _, _ in asm_rows if t["kind"] == kind}
    kp = {p for p in allp if alone.get((p, "0"), 0) >= KNOWN and alone.get((p, "100"), 0) >= KNOWN}
    cells = []
    for s in shares:
        v = [x for t, x, _ in asm_rows if t["kind"] == kind and t["pair"] in kp and float(t["b_pct"]) == s]
        cells.append(f"{sum(rule(x) for x in v)}/{len(v)}")
    print(f"| {kind} | {len(kp)} ({len(allp)}) | " + " | ".join(cells) + " |")
KNOWN_PAIRS = {p for p in {t["pair"] for t, _, _ in asm_rows} if alone.get((p, "0"), 0) >= KNOWN and alone.get((p, "100"), 0) >= KNOWN}

print("\n## Detection rate by minority share (read mixtures, default settings); parent B = minority share 100 - A share\n")
for depth in ("low", "high"):
    shs = sorted({100 - pa for _, d, pa, _, _ in reads_rows if d == depth and pa not in (100.0, 0.0)})
    print(f"\n{depth} depth (4x = low, 20x = high)\n")
    print("| pair class | n pairs | " + " | ".join(f"{s:g}%" for s in shs) + " |")
    print("|---|---|" + "---|" * len(shs))
    for dist in ("same_species", "same_genus", "same_order", "diff_phylum"):
        pairs = {r["id"].split("__")[0] for d2, d, pa, x, r in reads_rows if d2 == dist and d == depth}
        cells = []
        for s in shs:
            v = [x for d2, d, pa, x, r in reads_rows if d2 == dist and d == depth and abs(100 - pa - s) < 1e-6]
            cells.append(f"{sum(rule(x) for x in v)}/{len(v)}")
        print(f"| {dist} | {len(pairs)} | " + " | ".join(cells) + " |")

print(f"\n### the same, only read pairs whose two parents are each explained >= {KNOWN:g} by their genus alone (and not themselves flagged)\n")
palone = {}
for d2, d, pa, x, r in reads_rows:
    if x and pa in (100.0, 0.0):
        palone[(r["id"].split("__")[0], d, pa)] = x
for depth in ("low", "high"):
    shs = sorted({100 - pa for _, d, pa, _, _ in reads_rows if d == depth and pa not in (100.0, 0.0)})
    print(f"\n{depth} depth\n")
    print("| pair class | n pairs (of all) | " + " | ".join(f"{s:g}%" for s in shs) + " |")
    print("|---|---|" + "---|" * len(shs))
    for dist in ("same_species", "same_genus", "same_order", "diff_phylum"):
        allp = {r["id"].split("__")[0] for d2, d, pa, x, r in reads_rows if d2 == dist and d == depth}
        kp = set()
        for pid in allp:
            a_, b_ = palone.get((pid, depth, 100.0)), palone.get((pid, depth, 0.0))
            if a_ and b_ and a_["e"] >= KNOWN and b_["e"] >= KNOWN and not rule(a_) and not rule(b_):
                kp.add(pid)
        cells = []
        for s in shs:
            v = [x for d2, d, pa, x, r in reads_rows if d2 == dist and d == depth and r["id"].split("__")[0] in kp and abs(100 - pa - s) < 1e-6]
            cells.append(f"{sum(rule(x) for x in v)}/{len(v)}")
        print(f"| {dist} | {len(kp)} ({len(allp)}) | " + " | ".join(cells) + " |")

print("\n## Abundance error of the explained fractions\n")
print("Assemblies: true share of B = bases of B / total. Reads: true share of B = read-pair share, which is a hash share only when the parents have the same genome size and the "
      "minority is above the abundance cutoff (reads are shown for 20x and shares >= 5% only).\n")
print("| set | class | true B share | n | detected | median fraction of the second taxon | median error (abs, points) | median explained fraction of A (true share of A) |")
print("|---|---|---|---|---|---|---|---|")


def row(label, cls, q, xs):
    xs = [x for x in xs if x]
    det = [x for x in xs if rule(x)]
    if not det:
        print(f"| {label} | {cls} | {q:g}% | {len(xs)} | 0 | | | |")
        return
    mf = statistics.median(x["f"] for x in det)
    err = statistics.median(abs(x["f"] - q / 100) for x in det) * 100
    ea = statistics.median(x["e"] for x in det)
    print(f"| {label} | {cls} | {q:g}% | {len(xs)} | {len(det)} | {mf:.4f} | {err:.2f} | {ea:.3f} ({1 - q / 100:.3f}) |")


for kind in ("same_order", "diff_phylum"):
    for s in (1, 2, 5, 10, 20, 30, 50):
        row("assemblies", kind, s, [x for t, x, _ in asm_rows if t["kind"] == kind and float(t["b_pct"]) == s])
for dist in ("same_order", "diff_phylum"):
    for s in (5, 10, 20, 50):
        row("reads 20x", dist, s, [x for d2, d, pa, x, r in reads_rows if d2 == dist and d == "high" and abs(100 - pa - s) < 1e-6])

if a.sweep:
    print("\n## Sweep: false-mixture rate (all negatives pooled) and detection (assemblies, different-genus classes, B >= 5%; reads 20x, B >= 10%)\n")
    allneg = [x for v in neg.values() for _, x in v if x]
    pos_a = [x for t, x, _ in asm_rows if t["kind"] != "same_genus" and 5 <= float(t["b_pct"]) < 100]
    pos_r = [x for d2, d, pa, x, r in reads_rows if d2 in ("same_order", "diff_phylum") and d == "high" and 100 - pa >= 10 and pa > 0]
    print("| min fraction | min hashes | min concentration | min of rest | false rate (all negatives) | asm detect | reads detect |")
    print("|---|---|---|---|---|---|---|")
    for f in (0.005, 0.01, 0.02):
        for h in (10, 20, 50):
            for c in (0.0, 0.5):
                for rs_ in (0.0, 0.2, 0.3, 0.5):
                    fr = sum(rule(x, f, h, c, rs_) for x in allneg) / max(1, len(allneg))
                    da = sum(rule(x, f, h, c, rs_) for x in pos_a if x) / max(1, len(pos_a))
                    dr = sum(rule(x, f, h, c, rs_) for x in pos_r if x) / max(1, len(pos_r))
                    print(f"| {f} | {h} | {c} | {rs_} | {100 * fr:.2f}% | {100 * da:.0f}% | {100 * dr:.0f}% |")
