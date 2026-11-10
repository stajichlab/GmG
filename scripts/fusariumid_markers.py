#!/usr/bin/env python3
"""GmG marker mode for FusariumID: build inputs for a marker index, and score held-out marker queries.

A marker index is an ordinary gmg index built at scaled=1 (every k-mer), with one reference "genome" per FusariumID
strain holding all of that strain's marker sequences (tef1, rpb1, rpb2, CaM, tub2, ITS, ...). A query is any subset of
markers (e.g. one Sanger tef1 read); gmg ranks strains by shared k-mers and estimates identity from the query's
containment, exactly as for genomes.

prepare RECORDS LOCI_DIR OUT [--holdout 0.2] [--seed 1] [--loci tef1,rpb1,rpb2,CaM,tub2,ITS]
    records.tsv and loci/*.fa from bench/F1_fusariumid/parse_fusariumid.py.
    Writes OUT/refs/<id>.fa (one per reference strain), OUT/meta.tsv (gmg build metadata; lineage
    Ascomycota/Sordariomycetes/Hypocreales/Nectriaceae/<genus>/<species>, label = strain), and with --holdout > 0,
    strains drawn at random (per seed) as queries: OUT/queries/<locus>/<seq_id>.fa plus OUT/queries_truth.tsv.
score RESULTS_TSV TRUTH_TSV COMPLEX_TSV OUT_TSV [--min-precision 0.95]
    RESULTS_TSV from `gmg run` on the queries with permissive thresholds (species min_ani 0). Sweeps the species ANI
    cutoff per locus and writes called/precision for species calls (answerable queries: species present among the
    reference strains) and the false-call rate on queries whose species is absent; complex calls (FOSC) are scored
    separately. Prints the lowest cutoff whose pooled precision is >= --min-precision.

The FusariumID data are curated by fusarium.org; do not publish an index built from them without their agreement.
"""
import argparse
import collections
import os
import random
import re
import sys

LINEAGE = ["Ascomycota", "Sordariomycetes", "Hypocreales", "Nectriaceae"]
CUTOFFS = [0.98, 0.985, 0.99, 0.992, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999, 0.9995, 1.0]


def read_tsv(path):
    with open(path) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(cols, line.rstrip("\n").split("\t")))


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def safe(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")


def prepare(a):
    loci = a.loci.split(",")
    rec = {r["seq_id"]: r for r in read_tsv(a.records)}
    seqs = {}
    for l in loci:
        p = os.path.join(a.loci_dir, f"{l}.fa")
        if os.path.exists(p):
            seqs.update(dict(read_fasta(p)))
    strains = collections.defaultdict(list)
    for sid in seqs:
        r = rec[sid]
        strains[(r["strain"], r["species"])].append(sid)
    keys = sorted(strains)
    rng = random.Random(a.seed)
    held = set(rng.sample(keys, int(round(a.holdout * len(keys))))) if a.holdout > 0 else set()
    os.makedirs(os.path.join(a.out, "refs"), exist_ok=True)
    with open(os.path.join(a.out, "meta.tsv"), "w") as meta:
        meta.write("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\tlabel\n")
        for key in keys:
            if key in held:
                continue
            strain, sp = key
            ids = sorted(strains[key], key=lambda s: int(s) if s.isdigit() else s)
            gid = f"fid{ids[0]}_{safe(strain)}"
            with open(os.path.join(a.out, "refs", f"{gid}.fa"), "w") as f:
                for sid in ids:
                    f.write(f">{sid} {rec[sid]['locus']}\n{seqs[sid]}\n")
            meta.write("\t".join([gid, "0", *LINEAGE, sp.split()[0], sp, f"{sp} {strain}"]) + "\n")
    if held:
        with open(os.path.join(a.out, "queries_truth.tsv"), "w") as t:
            t.write("query_id\tlocus\tspecies\tstrain\tcomplex\n")
            for key in sorted(held):
                strain, sp = key
                for sid in strains[key]:
                    l = rec[sid]["locus"]
                    d = os.path.join(a.out, "queries", l)
                    os.makedirs(d, exist_ok=True)
                    qid = f"q{sid}"
                    with open(os.path.join(d, f"{qid}.fa"), "w") as f:
                        f.write(f">{qid} {l}\n{seqs[sid]}\n")
                    t.write(f"{qid}\t{l}\t{sp}\t{strain}\t{rec[sid]['complex']}\n")
    print(f"{len(keys) - len(held)} reference strains, {len(held)} held-out strains, loci {loci}", file=sys.stderr)


def load_complexes(path):
    t = {}
    for line in open(path):
        if line.startswith("#") or line.startswith("species\t"):
            continue
        sp, cx, col = line.rstrip("\n").split("\t")
        t[sp] = (cx, col == "1")
    return t


def score(a):
    truth = {r["query_id"]: r for r in read_tsv(a.truth)}
    cx = load_complexes(a.complexes)
    res = [r for r in read_tsv(a.results) if r["query_id"] in truth]
    ref_species = {r["species"] for r in read_tsv(a.meta)}
    rows = []
    groups = collections.defaultdict(list)
    for r in res:
        t = truth[r["query_id"]]
        groups[t["locus"]].append((r, t))
        groups["all"].append((r, t))
    best = None
    for locus, items in sorted(groups.items(), key=lambda kv: (kv[0] == "all", kv[0])):
        g_called = sum(bool(r["call_genus"]) for r, _ in items if "call_genus" in r)
        g_right = sum(r["call_genus"] == t["species"].split()[0] for r, t in items if r.get("call_genus"))
        for c in CUTOFFS:
            n_ans = n_abs = called = right = false_abs = cx_called = cx_right = 0
            for r, t in items:
                ani = float(r["ani"] or 0)
                answerable = t["species"] in ref_species
                n_ans += answerable
                n_abs += not answerable
                complex_ = r.get("species_complex", "")
                collapsed = bool(complex_) and "species_complex" in r.get("withheld", "")
                if collapsed:
                    cx_called += 1
                    cx_right += cx.get(t["species"], ("", False))[0] == complex_
                    continue
                sp = r["call_species"]
                if sp and ani >= c:
                    if answerable:
                        called += 1
                        right += sp == t["species"]
                    else:
                        false_abs += 1
            row = dict(
                locus=locus, cutoff=c, queries=len(items), answerable=n_ans, absent=n_abs,
                species_called=f"{called / max(n_ans, 1):.4f}", species_precision=f"{right / max(called, 1):.4f}",
                absent_false_call=f"{false_abs / max(n_abs, 1):.4f}",
                complex_calls=cx_called, complex_precision=f"{cx_right / max(cx_called, 1):.4f}",
                genus_called=f"{g_called / max(len(items), 1):.4f}", genus_precision=f"{g_right / max(g_called, 1):.4f}",
            )
            rows.append(row)
            if locus == "all" and best is None and called and right / called >= a.min_precision:
                best = c
    with open(a.out, "w") as f:
        cols = list(rows[0].keys())
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"lowest species ANI cutoff with pooled precision >= {a.min_precision}: {best}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("records")
    p.add_argument("loci_dir")
    p.add_argument("out")
    p.add_argument("--holdout", type=float, default=0.0)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--loci", default="tef1,rpb1,rpb2,CaM,tub2,ITS")
    s = sub.add_parser("score")
    s.add_argument("results")
    s.add_argument("truth")
    s.add_argument("complexes")
    s.add_argument("out")
    s.add_argument("--meta", required=True, help="meta.tsv of the reference strains (defines which species are answerable)")
    s.add_argument("--min-precision", type=float, default=0.95)
    a = ap.parse_args()
    prepare(a) if a.cmd == "prepare" else score(a)


if __name__ == "__main__":
    main()
