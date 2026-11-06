#!/usr/bin/env python3
"""Compare the sketch-level answer from subsamples with the whole run (output of 06_subsample_test.sbatch).
usage (from GmG_popgen): 07_subsample_compare.py subtest SLUG [SLUG ...]   (BIN=path to popsketch)
For each variant (full, stride1g, stride300, head300, head100) popsketch is run on the folder of sketches with --flagged mark (so all samples
get distances), then against `full`: distances (1 - ANI) of all sample pairs (Pearson r, mean and max absolute difference, and the mean
difference as a multiple of the variant's noise floor), nearest neighbour agreement, and the coverage and size estimates per sample."""
import csv, json, math, os, subprocess, sys

BIN = os.environ.get("BIN", "popsketch_build/target/release/popsketch")
VARIANTS = os.environ.get("VARIANTS", "full,stride1g,stride300,head300,head100").split(",")
# --from-sketches wants the parameters the sketches were made with
MAX_BASES = {"full": 0, "stride1g": 1_000_000_000, "stride300": 300_000_000, "head300": 300_000_000, "head100": 100_000_000}


def run(base, slug, v):
    out = f"{base}/{slug}/_out_{v}"
    subprocess.run([BIN, "run", "--from-sketches", f"{base}/{slug}/{v}", "-o", out, "--max-bases", str(MAX_BASES.get(v, 0)), "--flagged", "mark", "--min-cov", "0", "--no-pca", "--no-tree", "--threads", "4", "-q"], check=True)
    return out


def ani(out):
    rows = list(csv.reader(open(f"{out}/ani.tsv"), delimiter="\t"))
    ids = rows[0][1:]
    return ids, {(a, b): 1 - float(x) for a, r in zip(ids, rows[1:]) for b, x in zip(ids, r[1:]) if a < b}


def qc(out):
    return {r["id"]: r for r in csv.DictReader(open(f"{out}/qc.tsv"), delimiter="\t")}


def pearson(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sx = math.sqrt(sum((a - mx) ** 2 for a in x)); sy = math.sqrt(sum((b - my) ** 2 for b in y))
    return sxy / (sx * sy) if sx and sy else float("nan")


def nn(ids, d):
    return {i: min((j for j in ids if j != i), key=lambda j: d[(min(i, j), max(i, j))]) for i in ids}


base = sys.argv[1]
for slug in sys.argv[2:]:
    outs = {v: run(base, slug, v) for v in VARIANTS}
    data = {v: ani(outs[v]) for v in VARIANTS}
    q = {v: qc(outs[v]) for v in VARIANTS}
    floor = {v: json.load(open(f"{outs[v]}/popsketch.json"))["noise_floor"]["value"] for v in VARIANTS}
    ids0, dfull = data["full"]
    ids = [i for i in ids0 if all(i in data[v][0] for v in VARIANTS)]
    pairs = [(a, b) for a in ids for b in ids if a < b]
    nfull = nn(ids, dfull)
    print(f"\n{slug}: {len(ids)} samples, {len(pairs)} pairs; noise floor (1 - ANI) full {floor['full']:.4f}")
    print(f"{'variant':10} {'floor':>7} {'r vs full':>9} {'mean |d|':>9} {'max |d|':>8} {'mean|d|/floor':>13} {'NN same':>8} {'cov med':>8} {'size med Mb':>11}")
    for v in VARIANTS:
        dv = data[v][1]
        x = [dfull[p] for p in pairs]; y = [dv[p] for p in pairs]
        diff = [abs(a - b) for a, b in zip(x, y)]
        nv = nn(ids, dv)
        same = sum(nv[i] == nfull[i] for i in ids) / len(ids)
        cov = sorted(float(q[v][i]["kmer_coverage"] or 0) for i in ids)[len(ids) // 2]
        size = sorted(float(q[v][i]["est_size_mb"] or 0) for i in ids)[len(ids) // 2]
        print(f"{v:10} {floor[v]:7.4f} {pearson(x, y):9.4f} {sum(diff)/len(diff):9.5f} {max(diff):8.5f} {sum(diff)/len(diff)/floor[v]:13.2f} {same:8.2f} {cov:8.1f} {size:11.1f}")
    # per-sample coverage and size against the whole run
    for v in VARIANTS[1:]:
        dc = [float(q[v][i]["kmer_coverage"] or 0) / max(float(q["full"][i]["kmer_coverage"] or 0), 1e-9) for i in ids]
        ds = [float(q[v][i]["est_size_mb"] or 0) / max(float(q["full"][i]["est_size_mb"] or 0), 1e-9) for i in ids]
        print(f"  {v:10} coverage / full: median {sorted(dc)[len(dc)//2]:.2f} range {min(dc):.2f}-{max(dc):.2f};  size / full: median {sorted(ds)[len(ds)//2]:.2f} range {min(ds):.2f}-{max(ds):.2f}")
