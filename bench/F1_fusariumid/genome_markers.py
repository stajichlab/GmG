#!/usr/bin/env python3
"""In-silico barcoding: call a genome's FusariumID species from its marker loci.

Input: blastn of FusariumID sequences (queries) against ONE genome (subject), -outfmt
    "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send bitscore"
plus FusariumID records.tsv. For each FusariumID sequence the best HSP (bitscore) is kept when it covers >= --min-cov of
the sequence; distance = (mismatches + gap openings) / aligned length.

Calls
    per locus: the species of the closest FusariumID sequence(s) (ties -> several species, ';'-joined), its distance and
               differences, and the distance of the closest sequence of any OTHER species (the margin).
    combined:  per FusariumID strain, sum(differences) / sum(aligned length) over the loci in --loci that hit; strains
               need >= 2 such loci; the closest strain's species is the call (ties -> several species).

usage: genome_markers.py GENOME_ID BLAST_TSV RECORDS_TSV [--loci tef1,rpb1,rpb2,CaM,tub2] [--header]
prints one TSV row per genome (and the header with --header)
"""
import argparse
import collections

TOL = 1e-12
COLS = ["genome", "n_loci_hit", "combined_species", "combined_diffs", "combined_len", "combined_dist",
        "combined_other_species", "combined_other_dist", "combined_margin"]


def read_records(path):
    rec = {}
    with open(path) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            rec[r["seq_id"]] = r
    return rec


def best_hits(path, min_cov):
    best = {}
    with open(path) as f:
        for line in f:
            q, _s, _p, ln, mm, go, _qs, _qe, ql, _ss, _se, bits = line.rstrip("\n").split("\t")[:12]
            ln, ql, b = int(ln), int(ql), float(bits)
            if ln < min_cov * ql:
                continue
            if q not in best or b > best[q][2]:
                best[q] = (int(mm) + int(go), ln, b)
    return best


def closest(items):
    """items: (species, dist, diffs, length). -> (species set, dist, diffs, length, other_species, other_dist)."""
    items = sorted(items, key=lambda x: (x[1], -x[3]))
    if not items:
        return set(), None, None, None, "", None
    d0 = items[0][1]
    top = {sp for sp, d, *_ in items if abs(d - d0) <= TOL}
    other = next(((sp, d) for sp, d, *_ in items if sp not in top), ("", None))
    return top, d0, items[0][2], items[0][3], other[0], other[1]


def fmt(x, nd=5):
    return "" if x is None else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def call(genome, blast_tsv, records, loci, min_cov=0.6):
    rec = records
    best = best_hits(blast_tsv, min_cov)
    per_locus = collections.defaultdict(list)
    per_strain = collections.defaultdict(dict)
    for q, (d, ln, _b) in best.items():
        r = rec.get(q)
        if r is None or r["locus"] not in loci:
            continue
        per_locus[r["locus"]].append((r["species"], d / ln, d, ln))
        key = (r["strain"], r["species"])
        if r["locus"] not in per_strain[key]:
            per_strain[key][r["locus"]] = (d, ln)
    row = {"genome": genome, "n_loci_hit": sum(1 for l in loci if per_locus.get(l))}
    comb = []
    for (strain, sp), m in per_strain.items():
        if len(m) >= 2:
            d = sum(v[0] for v in m.values())
            ln = sum(v[1] for v in m.values())
            comb.append((sp, d / ln, d, ln))
    top, dist, diffs, ln, osp, od = closest(comb)
    row.update(
        combined_species=";".join(sorted(top)),
        combined_diffs=fmt(diffs),
        combined_len=fmt(ln),
        combined_dist=fmt(dist),
        combined_other_species=osp,
        combined_other_dist=fmt(od),
        combined_margin=fmt(None if od is None or dist is None else od - dist),
    )
    for l in loci:
        top, dist, diffs, ln, osp, od = closest(per_locus.get(l, []))
        row[f"{l}_species"] = ";".join(sorted(top))
        row[f"{l}_dist"] = fmt(dist)
        row[f"{l}_diffs"] = fmt(diffs)
        row[f"{l}_other_dist"] = fmt(od)
    return row


def header(loci):
    return COLS + [f"{l}_{c}" for l in loci for c in ("species", "dist", "diffs", "other_dist")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome")
    ap.add_argument("blast")
    ap.add_argument("records")
    ap.add_argument("--loci", default="tef1,rpb1,rpb2,CaM,tub2")
    ap.add_argument("--min-cov", type=float, default=0.6)
    ap.add_argument("--header", action="store_true")
    a = ap.parse_args()
    loci = a.loci.split(",")
    cols = header(loci)
    if a.header:
        print("\t".join(cols))
    row = call(a.genome, a.blast, read_records(a.records), loci, a.min_cov)
    print("\t".join(str(row.get(c, "")) for c in cols))


if __name__ == "__main__":
    main()
