#!/usr/bin/env python3
"""How much marker sequence separates fusarioid species? Leave-one-out on the FusariumID database.

Inputs (from parse_fusariumid.py and run_markers.sbatch):
    FID/records.tsv               one row per sequence
    FID/loci/<locus>.fa           sequences per locus
    BLAST/<locus>.tsv             all-vs-all blastn, -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen slen bitscore"

Method
    distance between two sequences = (mismatches + gap openings) / aligned length of their best HSP, kept when the HSP
    covers >= 200 bp and >= 50% of the shorter sequence. A query is compared only with *other strains*. The call is the
    species of the nearest neighbour(s): correct (one species, the right one), ambiguous (tie between the right and
    another species), wrong (the right species is not among the nearest), no hit. A query is answerable when its species
    has >= 2 strains with data for the loci in question; singletons are reported apart (they can only be "novel").

    Multilocus: per strain pair, distance = sum(differences) / sum(aligned length) over the loci of the set; both strains
    must have every locus of the set (complete-data subset, so the n shrinks as loci are added).

    k-mer variant (GmG's currency): canonical k-mers (k=21), FracMinHash at scaled s (hash < 2^64/s), Jaccard, nearest
    neighbour as above. At s=1 this is all k-mers; GmG's genome sketches use s >= 1000.

Outputs (OUTDIR): loci.tsv, multilocus.tsv, complex.tsv, nearest_species.tsv, kmer.tsv, per_query.tsv
usage: marker_resolution.py FID BLAST OUTDIR [--kmer-loci tef1,rpb2,...] [--scaled 1,10,30,100]
"""
import argparse
import collections
import os
import sys

MIN_ALN = 200
MIN_COV = 0.5
TOL = 1e-12


def read_records(path):
    rec = {}
    with open(path) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            r["strain_key"] = (r["strain"], r["species"])
            rec[r["seq_id"]] = r
    return rec


def read_blast(path):
    """Best HSP per (query, subject) pair -> {q: {s: (differences, aligned_length)}}."""
    best = {}
    with open(path) as f:
        for line in f:
            q, s, _pid, ln, mm, go, ql, sl, bits = line.split("\t")[:9]
            if q == s:
                continue
            ln, ql, sl = int(ln), int(ql), int(sl)
            if ln < MIN_ALN or ln < MIN_COV * min(ql, sl):
                continue
            b = float(bits)
            d = best.setdefault(q, {})
            if s not in d or b > d[s][2]:
                d[s] = (int(mm) + int(go), ln, b)
    return {q: {s: (v[0], v[1]) for s, v in d.items()} for q, d in best.items()}


def outcome(truth, nn_species):
    if not nn_species:
        return "nohit"
    if nn_species == {truth}:
        return "correct"
    return "ambiguous" if truth in nn_species else "wrong"


def nearest(items, truth):
    """items: iterable of (species, distance, differences). Returns (nn_species set, nn_dist, nn_con, nn_het, nn_het_diff)."""
    nn, nn_sp = float("inf"), set()
    con, het, het_d = float("inf"), float("inf"), None
    for sp, p, diff in items:
        if p < nn - TOL:
            nn, nn_sp = p, {sp}
        elif abs(p - nn) <= TOL:
            nn_sp.add(sp)
        if sp == truth:
            con = min(con, p)
        elif p < het or (p == het and (het_d is None or diff < het_d)):
            het, het_d = p, diff
    return nn_sp, nn, con, het, het_d


def summarise(rows, key_cols):
    """rows: dicts with outcome, answerable, gap. Returns summary dict."""
    n = len(rows)
    ans = [r for r in rows if r["answerable"]]
    c = collections.Counter(r["outcome"] for r in ans)
    hd = sorted(r["het_diff"] for r in rows if r["het_diff"] is not None)
    out = dict(key_cols)
    out.update(
        n=n,
        answerable=len(ans),
        correct=c["correct"],
        ambiguous=c["ambiguous"],
        wrong=c["wrong"],
        nohit=c["nohit"],
        accuracy=f"{c['correct'] / len(ans):.4f}" if ans else "",
        barcode_gap=f"{sum(r['gap'] for r in ans) / len(ans):.4f}" if ans else "",
        median_diffs_to_other_species=hd[len(hd) // 2] if hd else "",
        frac_le2_diffs=f"{sum(x <= 2 for x in hd) / len(hd):.4f}" if hd else "",
    )
    return out


def write_tsv(path, rows):
    if not rows:
        return
    cols = list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def per_locus(rec, hits, locus):
    seqs = [s for s, r in rec.items() if r["locus"] == locus]
    strains_per_sp = collections.Counter()
    for sp, _st in {(rec[s]["species"], rec[s]["strain_key"]) for s in seqs}:
        strains_per_sp[sp] += 1
    rows = []
    for q in seqs:
        r = rec[q]
        truth = r["species"]
        items = [
            (rec[s]["species"], d / ln, d)
            for s, (d, ln) in hits.get(q, {}).items()
            if rec[s]["strain_key"] != r["strain_key"]
        ]
        nn_sp, nn, con, het, het_d = nearest(items, truth)
        rows.append(
            dict(
                level="locus",
                set=locus,
                query=q,
                species=truth,
                complex=r["complex"],
                answerable=strains_per_sp[truth] >= 2,
                outcome=outcome(truth, nn_sp),
                nn_dist=nn,
                gap=con < het,
                het_diff=het_d,
            )
        )
    return rows


def strain_pairs(rec, hits_by_locus):
    """{locus: {strain_a: {strain_b: (d, L)}}} using the first sequence of a strain per locus."""
    first = {}
    for s, r in rec.items():
        first.setdefault((r["strain_key"], r["locus"]), s)
    keep = set(first.values())
    out = {}
    for locus, hits in hits_by_locus.items():
        m = {}
        for q, d in hits.items():
            if q not in keep:
                continue
            a = rec[q]["strain_key"]
            row = m.setdefault(a, {})
            for s, v in d.items():
                if s in keep:
                    row[rec[s]["strain_key"]] = v
        out[locus] = m
    return out


def multilocus(rec, sp_hits, combo):
    strains_loci = collections.defaultdict(set)
    meta = {}
    for r in rec.values():
        strains_loci[r["strain_key"]].add(r["locus"])
        meta[r["strain_key"]] = r
    full = {st for st, ls in strains_loci.items() if all(l in ls for l in combo)}
    per_sp = collections.Counter(st[1] for st in full)
    rows, nearest_other = [], {}
    for a in full:
        truth = a[1]
        cand = None
        for l in combo:
            ks = set(sp_hits[l].get(a, {}))
            cand = ks if cand is None else cand & ks
        items = []
        for b in cand or ():
            if b not in full or b == a:
                continue
            d = sum(sp_hits[l][a][b][0] for l in combo)
            ln = sum(sp_hits[l][a][b][1] for l in combo)
            items.append((b[1], d / ln, d))
        nn_sp, nn, con, het, het_d = nearest(items, truth)
        if het_d is not None:
            nearest_other[truth] = min(nearest_other.get(truth, het_d), het_d)
        rows.append(
            dict(
                level="multilocus",
                set="+".join(combo),
                query=f"{a[0]}",
                species=truth,
                complex=meta[a]["complex"],
                answerable=per_sp[truth] >= 2,
                outcome=outcome(truth, nn_sp),
                nn_dist=nn,
                gap=con < het,
                het_diff=het_d,
            )
        )
    return rows, nearest_other


# ---------- k-mer variant ----------


def kmer_hashes(seq, k):
    import numpy as np

    code = np.frombuffer(seq.encode(), dtype=np.uint8)
    lut = np.full(256, 4, dtype=np.uint8)
    for i, c in enumerate(b"ACGT"):
        lut[c] = i
    c = lut[code]
    if len(c) < k:
        return np.zeros(0, dtype=np.uint64)
    from numpy.lib.stride_tricks import sliding_window_view

    w = sliding_window_view(c, k)
    ok = (w < 4).all(axis=1)
    w = w[ok].astype(np.uint64)
    if len(w) == 0:
        return np.zeros(0, dtype=np.uint64)
    pw = (np.uint64(4) ** np.arange(k - 1, -1, -1, dtype=np.uint64)).astype(np.uint64)
    fw = (w * pw).sum(axis=1, dtype=np.uint64)
    rc = ((np.uint64(3) - w[:, ::-1]) * pw).sum(axis=1, dtype=np.uint64)
    can = np.minimum(fw, rc)
    # splitmix64 finaliser as the FracMinHash hash
    z = can + np.uint64(0x9E3779B97F4A7C15)
    z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    z = z ^ (z >> np.uint64(31))
    return np.unique(z)


def kmer_loo(rec, fasta, locus, k, scaled_list):
    import numpy as np
    from scipy import sparse

    ids, hs = [], []
    name = None
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].strip()
            else:
                ids.append(name)
                hs.append(kmer_hashes(line.strip(), k))
    strains_per_sp = collections.Counter({})
    for sp, _st in {(rec[s]["species"], rec[s]["strain_key"]) for s in ids}:
        strains_per_sp[sp] += 1
    out = []
    for s in scaled_list:
        thr = np.uint64((2**64 - 1) // s)
        sub = [h[h <= thr] for h in hs]
        sizes = np.array([len(h) for h in sub])
        allh = np.concatenate(sub) if sub else np.zeros(0, dtype=np.uint64)
        uniq, inv = np.unique(allh, return_inverse=True)
        rows = np.repeat(np.arange(len(sub)), sizes)
        X = sparse.csr_matrix((np.ones(len(inv), dtype=np.float32), (rows, inv)), shape=(len(sub), len(uniq)))
        S = (X @ X.T).tocsr()
        res = []
        for i, q in enumerate(ids):
            r = rec[q]
            truth = r["species"]
            lo, hi = S.indptr[i], S.indptr[i + 1]
            items = []
            for j, sh in zip(S.indices[lo:hi], S.data[lo:hi]):
                if j == i or rec[ids[j]]["strain_key"] == r["strain_key"]:
                    continue
                jac = sh / (sizes[i] + sizes[j] - sh)
                items.append((rec[ids[j]]["species"], 1.0 - jac, None))
            nn_sp, *_ = nearest([(a, b, 0) for a, b, _ in items], truth)
            res.append(dict(answerable=strains_per_sp[truth] >= 2, outcome=outcome(truth, nn_sp), gap=False, het_diff=None))
        sm = summarise(res, dict(locus=locus, k=k, scaled=s))
        sm["mean_hashes_per_seq"] = f"{sizes.mean():.1f}"
        for c in ("barcode_gap", "median_diffs_to_other_species", "frac_le2_diffs"):
            sm.pop(c, None)
        out.append(sm)
        print(f"kmer {locus} s={s}: acc={sm['accuracy']} hashes/seq={sm['mean_hashes_per_seq']}", file=sys.stderr)
    return out


COMBOS = [
    ["tef1"],
    ["rpb2"],
    ["rpb1"],
    ["CaM"],
    ["tub2"],
    ["ITS"],
    ["LSU"],
    ["tef1", "rpb2"],
    ["tef1", "rpb1"],
    ["rpb1", "rpb2"],
    ["tef1", "rpb1", "rpb2"],
    ["tef1", "rpb2", "CaM"],
    ["tef1", "rpb2", "CaM", "tub2"],
    ["tef1", "rpb1", "rpb2", "CaM", "tub2"],
]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fid")
    ap.add_argument("blast")
    ap.add_argument("outdir")
    ap.add_argument("--kmer-loci", default="tef1,rpb2,rpb1,CaM,tub2,ITS")
    ap.add_argument("--scaled", default="1,3,10,30,100")
    ap.add_argument("-k", type=int, default=21)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    rec = read_records(os.path.join(a.fid, "records.tsv"))
    loci = sorted({r["locus"] for r in rec.values()})
    hits = {}
    for l in loci:
        p = os.path.join(a.blast, f"{l}.tsv")
        if os.path.exists(p):
            hits[l] = read_blast(p)
            print(f"blast {l}: {len(hits[l])} queries with hits", file=sys.stderr)

    all_rows, loc_sum, cx_rows = [], [], []
    for l in hits:
        rows = per_locus(rec, hits[l], l)
        all_rows += rows
        loc_sum.append(summarise(rows, dict(set=l)))
    write_tsv(os.path.join(a.outdir, "loci.tsv"), loc_sum)

    sp_hits = strain_pairs(rec, hits)
    ml_sum, nearest_rows = [], []
    for combo in COMBOS:
        if not all(l in sp_hits for l in combo):
            continue
        rows, near = multilocus(rec, sp_hits, combo)
        all_rows += rows
        ml_sum.append(summarise(rows, dict(set="+".join(combo), loci=len(combo))))
        for sp, d in sorted(near.items()):
            nearest_rows.append(dict(set="+".join(combo), species=sp, min_diffs_to_other_species=d))
        if combo in (["tef1"], ["tef1", "rpb1", "rpb2"], ["tef1", "rpb2", "CaM", "tub2"]):
            by_cx = collections.defaultdict(list)
            for r in rows:
                by_cx[r["complex"]].append(r)
            for cx, rs in sorted(by_cx.items(), key=lambda x: -len(x[1])):
                cx_rows.append(summarise(rs, dict(set="+".join(combo), complex=cx)))
    write_tsv(os.path.join(a.outdir, "multilocus.tsv"), ml_sum)
    write_tsv(os.path.join(a.outdir, "complex.tsv"), cx_rows)
    write_tsv(os.path.join(a.outdir, "nearest_species.tsv"), nearest_rows)
    write_tsv(
        os.path.join(a.outdir, "per_query.tsv"),
        [{k: v for k, v in r.items() if k != "gap"} | {"gap": int(r["gap"])} for r in all_rows],
    )

    km = []
    scaled = [int(x) for x in a.scaled.split(",") if x]
    for l in [x for x in a.kmer_loci.split(",") if x]:
        fa = os.path.join(a.fid, "loci", f"{l}.fa")
        if os.path.exists(fa):
            km += kmer_loo(rec, fa, l, a.k, scaled)
    write_tsv(os.path.join(a.outdir, "kmer.tsv"), km)


if __name__ == "__main__":
    main()
