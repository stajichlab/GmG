#!/usr/bin/env python3
"""Effect of gmg's species-complex rule (schema 10) on the F1 fusarioid benchmark, without re-running gmg.

Applies the same rule as gmg-core (`call.rs`): a species call needs top ANI >= 0.9938 (the default species gate; other
gates ignored, as in genome_resolution.py); when the species vote is in a collapsing complex of
crates/gmg-core/data/species_complexes.tsv and the top ANI >= 0.98, the complex is reported instead.
Scored against the FusariumID labels (genome_results/labels.tsv): a species call is right when the top hit's FusariumID
label equals the query's; a complex call is right when the query's FusariumID species is in that complex.

Also checks COMPLEX_MIN_ANI: for genomes whose FusariumID species is NOT in a collapsing complex, the highest ANI to any
genome that IS (from ani_s1000.tsv) should stay below 0.98.

usage: complex_effect.py W COMPLEX_TSV [--scaled 10000]   prints TSV-ish lines
"""
import argparse
import collections
import os

SPECIES_ANI, COMPLEX_ANI = 0.9938, 0.98


def read_tsv(path):
    with open(path) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(cols, line.rstrip("\n").split("\t")))


def load_complexes(path):
    t = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.startswith("species\t"):
                continue
            sp, cx, col = line.rstrip("\n").split("\t")
            t[sp] = (cx, col == "1")
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("w")
    ap.add_argument("complexes")
    ap.add_argument("--scaled", type=int, default=10000)
    a = ap.parse_args()
    cx = load_complexes(a.complexes)
    fid = {r["genome"]: r["fid_label"] for r in read_tsv(os.path.join(a.w, "genome_results", "labels.tsv")) if r["fid_label"]}
    count = collections.Counter(fid.values())

    def in_collapse(label):  # a FusariumID label (possibly "A;B") is in a collapsing complex if its first name is
        e = cx.get(label.split(";")[0])
        return e[0] if e and e[1] else None

    before = collections.Counter()
    after = collections.Counter()
    for r in read_tsv(os.path.join(a.w, "bench", f"fus_loo_s{a.scaled}.tsv")):
        q, t = fid.get(r["id"]), fid.get(r["top_id"])
        if not q or not t or count[q] < 2:
            continue
        ani = float(r["ani"])
        e = cx.get(r["species_vote"])
        collapsed = bool(e and e[1] and ani >= COMPLEX_ANI and r["genus_vote"])
        sp_call = ani >= SPECIES_ANI
        before["n"] += 1
        after["n"] += 1
        if sp_call:
            before["species_called"] += 1
            before["species_right"] += q == t
        if collapsed:
            after["complex_called"] += 1
            after["complex_right"] += in_collapse(q) == e[0]
        elif sp_call:
            after["species_called"] += 1
            after["species_right"] += q == t
    for name, c in (("before", before), ("after", after)):
        line = f"{name}\tn={c['n']}\tspecies_called={c['species_called']}\tspecies_precision={c['species_right'] / max(c['species_called'], 1):.4f}"
        if name == "after":
            line += f"\tcomplex_called={c['complex_called']}\tcomplex_precision={c['complex_right'] / max(c['complex_called'], 1):.4f}"
        print(line)

    # COMPLEX_MIN_ANI check
    worst = collections.defaultdict(float)
    best_in = collections.defaultdict(float)
    for r in read_tsv(os.path.join(a.w, "ani_s1000.tsv")):
        la, lb = fid.get(r["a"]), fid.get(r["b"])
        if not la or not lb:
            continue
        ani = float(r["ani"])
        ca, cb = in_collapse(la), in_collapse(lb)
        for g, mine, other in ((r["a"], ca, cb), (r["b"], cb, ca)):
            if mine is None and other is not None:
                worst[g] = max(worst[g], ani)
            if mine is not None and other is not None:
                best_in[g] = max(best_in[g], ani)
    out = sorted(worst.values(), reverse=True)
    n_out = sum(1 for g, l in fid.items() if in_collapse(l) is None)
    print(f"non-collapse genomes: {n_out}; with any ANI to a collapse-complex genome: {len(out)}; max {out[0] if out else 0:.4f};"
          f" >= {COMPLEX_ANI}: {sum(x >= COMPLEX_ANI for x in out)}")
    ins = sorted(best_in.values())
    print(f"collapse-complex genomes: {len(ins)}; nearest same-complex ANI min {ins[0] if ins else 0:.4f};"
          f" < {COMPLEX_ANI}: {sum(x < COMPLEX_ANI for x in ins)}")


if __name__ == "__main__":
    main()
