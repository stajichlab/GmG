#!/usr/bin/env python3
"""Join the F1 genome results: FusariumID marker labels vs NCBI labels, gmg leave-out accuracy, genome ANI gaps.

Inputs (build dir W):
    W/genomes/meta.tsv          index metadata of the fusarioid genomes (NCBI labels)
    W/bridge/calls_*.tsv        genome_markers.py rows (FusariumID call per genome)
    W/ani_s1000.tsv             genome_ani.py all pairs
    W/bench/{all,fus}_{loo,los}_s*.tsv   gmg bench outputs

FusariumID label of a genome ("fid"): the combined multilocus call when it is within --max-dist (default 0.01, i.e.
>= 99% identity over the concatenated loci); a tie between species keeps all names ("A;B"), so tied genomes form their
own group. Genomes with no call or a more distant one have no fid label and are not scored on it.

Outputs (OUTDIR):
    labels.tsv          per genome: NCBI species, fid label, distance, margin, agreement class
    label_summary.tsv   counts of agreement classes; top NCBI -> fid relabelings
    bench_summary.tsv   per reference set x mode x scaled: species accuracy vs NCBI and vs fid labels
    ani_gap.tsv         per labelling (ncbi, fid): barcode-gap fraction and quantiles of the nearest other-species ANI
    ani_species_pairs.tsv  closest pairs of distinct fid species: max genome ANI
usage: genome_resolution.py W OUTDIR
"""
import argparse
import collections
import glob
import os
import re

SPECIES_ANI = 0.9938  # gmg_core::call::Thresholds::default().min_ani[5]


def read_tsv(path):
    with open(path) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(cols, line.rstrip("\n").split("\t")))


def write_tsv(path, rows, cols=None):
    if not rows:
        return
    cols = cols or list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def q(xs, p):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(p * len(xs)))] if xs else ""


def fid_labels(rows, max_dist):
    out = {}
    for r in rows:
        sp, d = r["combined_species"], r["combined_dist"]
        if sp and d and float(d) <= max_dist:
            out[r["genome"]] = sp
    return out


def agreement(ncbi, fid):
    if not fid:
        return "no_fid_call"
    names = fid.split(";")
    if ncbi in names:
        return "same" if len(names) == 1 else "tie_includes_ncbi"
    if ncbi.split()[-1] in ("sp.",) or " cf. " in ncbi or ncbi.endswith(" sp"):
        return "ncbi_unnamed"
    return "differs"


def ani_gap(pairs, labels):
    """Per genome: max ANI to same label and to a different label (labels: id -> label; unlabelled ids ignored)."""
    con = collections.defaultdict(float)
    het = collections.defaultdict(float)
    for a, b, ani in pairs:
        la, lb = labels.get(a), labels.get(b)
        if la is None or lb is None:
            continue
        d = con if la == lb else het
        d[a] = max(d[a], ani)
        d[b] = max(d[b], ani)
    count = collections.Counter(labels.values())
    rows = []
    for g, l in labels.items():
        if count[l] < 2:
            continue
        rows.append((g, l, con.get(g, 0.0), het.get(g, 0.0)))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("w")
    ap.add_argument("outdir")
    ap.add_argument("--max-dist", type=float, default=0.01)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    meta = {r["id"]: r for r in read_tsv(os.path.join(a.w, "genomes", "meta.tsv"))}
    ncbi = {g: r["species"] for g, r in meta.items()}
    calls = []
    for p in sorted(glob.glob(os.path.join(a.w, "bridge", "calls_*.tsv"))):
        calls += list(read_tsv(p))
    calls_by = {r["genome"]: r for r in calls}
    fid = fid_labels(calls, a.max_dist)

    # 1. labels
    lab_rows, cls = [], collections.Counter()
    relabel = collections.Counter()
    for g in meta:
        c = calls_by.get(g, {})
        cl = agreement(ncbi[g], fid.get(g, ""))
        cls[cl] += 1
        if cl in ("differs", "ncbi_unnamed"):
            relabel[(ncbi[g], fid[g])] += 1
        lab_rows.append(
            dict(
                genome=g,
                ncbi_species=ncbi[g],
                fid_label=fid.get(g, ""),
                fid_call=c.get("combined_species", ""),
                fid_dist=c.get("combined_dist", ""),
                fid_diffs=c.get("combined_diffs", ""),
                fid_len=c.get("combined_len", ""),
                fid_margin=c.get("combined_margin", ""),
                fid_other=c.get("combined_other_species", ""),
                n_loci_hit=c.get("n_loci_hit", ""),
                agreement=cl,
            )
        )
    write_tsv(os.path.join(a.outdir, "labels.tsv"), lab_rows)
    summ = [dict(item="agreement", key=k, n=v) for k, v in cls.most_common()]
    tie = sum(1 for v in fid.values() if ";" in v)
    summ.append(dict(item="fid_labelled", key="all", n=len(fid)))
    summ.append(dict(item="fid_labelled", key="tie", n=tie))
    summ.append(dict(item="distinct_labels", key="ncbi", n=len(set(ncbi.values()))))
    summ.append(dict(item="distinct_labels", key="fid", n=len(set(fid.values()))))
    for (n_, f_), v in relabel.most_common(40):
        summ.append(dict(item="relabel", key=f"{n_} -> {f_}", n=v))
    write_tsv(os.path.join(a.outdir, "label_summary.tsv"), summ)

    # 2. bench
    bench_rows = []
    fid_count = collections.Counter(fid.values())
    for p in sorted(glob.glob(os.path.join(a.w, "bench", "*_s*.tsv"))):
        m = re.match(r"(all|fus)_(loo|los)_s(\d+)\.tsv", os.path.basename(p))
        if not m:
            continue
        refset, mode, scaled = m.group(1), m.group(2), int(m.group(3))
        rows = list(read_tsv(p))
        n = len(rows)
        ncbi_count = collections.Counter(ncbi[r["id"]] for r in rows)
        if mode == "loo":
            ans = [r for r in rows if ncbi_count[r["species_true"]] >= 2]
            top_ok = sum(r["species_top"] == r["species_true"] for r in ans)
            vote_ok = sum(r["species_vote"] == r["species_true"] for r in ans)
            called = sum(bool(r["species_vote"]) for r in ans)
            fa = [r for r in rows if fid.get(r["id"]) and fid_count[fid[r["id"]]] >= 2 and fid.get(r["top_id"])]
            fid_ok = sum(fid[r["top_id"]] == fid[r["id"]] for r in fa)
            # bench writes the raw vote; gmg's species call also needs top ANI >= SPECIES_ANI (other gates ignored here)
            g_ans = [r for r in ans if float(r["ani"]) >= SPECIES_ANI]
            g_ok = sum(r["species_top"] == r["species_true"] for r in g_ans)
            g_fa = [r for r in fa if float(r["ani"]) >= SPECIES_ANI]
            g_fok = sum(fid[r["top_id"]] == fid[r["id"]] for r in g_fa)
            bench_rows.append(
                dict(
                    refset=refset, mode=mode, scaled=scaled, n=n,
                    median_query_hashes=q([int(r["query_hashes"]) for r in rows], 0.5),
                    ncbi_answerable=len(ans),
                    ncbi_top_correct=f"{top_ok / len(ans):.4f}" if ans else "",
                    ncbi_vote_correct=f"{vote_ok / len(ans):.4f}" if ans else "",
                    ncbi_vote_called=f"{called / len(ans):.4f}" if ans else "",
                    fid_scored=len(fa),
                    fid_top_correct=f"{fid_ok / len(fa):.4f}" if fa else "",
                    gated_called=f"{len(g_ans) / len(ans):.4f}" if ans else "",
                    gated_ncbi_precision=f"{g_ok / len(g_ans):.4f}" if g_ans else "",
                    gated_fid_precision=f"{g_fok / len(g_fa):.4f}" if g_fa else "",
                    median_top_ani=q([float(r["ani"]) for r in rows], 0.5),
                )
            )
        else:  # los: the true species is gone; passing the species ANI gate would be a false species call
            false_sp = sum(float(r["ani"]) >= SPECIES_ANI for r in rows)
            genus_ok = sum(r["genus_vote"] == r["genus_true"] for r in rows)
            bench_rows.append(
                dict(
                    refset=refset, mode=mode, scaled=scaled, n=n,
                    median_query_hashes=q([int(r["query_hashes"]) for r in rows], 0.5),
                    los_false_species_calls=f"{false_sp / n:.4f}" if n else "",
                    los_genus_correct=f"{genus_ok / n:.4f}" if n else "",
                    median_top_ani=q([float(r["ani"]) for r in rows], 0.5),
                    los_top_ani_q90=q([float(r["ani"]) for r in rows], 0.9),
                )
            )
    cols = ["refset", "mode", "scaled", "n", "median_query_hashes", "ncbi_answerable", "ncbi_top_correct",
            "ncbi_vote_correct", "ncbi_vote_called", "fid_scored", "fid_top_correct", "gated_called",
            "gated_ncbi_precision", "gated_fid_precision", "los_false_species_calls",
            "los_genus_correct", "median_top_ani", "los_top_ani_q90"]
    write_tsv(os.path.join(a.outdir, "bench_summary.tsv"), bench_rows, cols)

    # 2b. per species complex (fid labels), fusarioid-only reference set
    cxmap = {}
    recp = os.path.join(a.w, "fid", "records.tsv")
    if os.path.exists(recp):
        for r in read_tsv(recp):
            cxmap.setdefault(r["species"], r["complex"].split(" (")[0])
    cx_rows = []
    for p in sorted(glob.glob(os.path.join(a.w, "bench", "fus_loo_s*.tsv"))):
        scaled = int(re.search(r"_s(\d+)\.tsv", p).group(1))
        lp = p.replace("_loo_", "_los_")
        los = {r["id"]: float(r["ani"]) for r in read_tsv(lp)} if os.path.exists(lp) else {}
        st = collections.defaultdict(lambda: [0, 0, 0, 0, 0])
        for r in read_tsv(p):
            f, t = fid.get(r["id"]), fid.get(r["top_id"])
            if not f or not t or fid_count[f] < 2:
                continue
            s = st[cxmap.get(f.split(";")[0], "?")]
            ok, gate = f == t, float(r["ani"]) >= SPECIES_ANI
            s[0] += 1
            s[1] += ok
            s[2] += gate
            s[3] += gate and ok
            s[4] += los.get(r["id"], 0.0) >= SPECIES_ANI
        for c, s in sorted(st.items(), key=lambda kv: -kv[1][0]):
            cx_rows.append(
                dict(scaled=scaled, complex=c, n=s[0], top_correct=f"{s[1] / s[0]:.4f}", gated_called=f"{s[2] / s[0]:.4f}",
                     gated_precision=f"{s[3] / s[2]:.4f}" if s[2] else "", los_above_gate=f"{s[4] / s[0]:.4f}")
            )
    write_tsv(os.path.join(a.outdir, "complex_bench.tsv"), cx_rows)

    # 3. ANI gap
    pairs = []
    ap_ = os.path.join(a.w, "ani_s1000.tsv")
    if os.path.exists(ap_):
        for r in read_tsv(ap_):
            pairs.append((r["a"], r["b"], float(r["ani"])))
    gap_rows, sp_pairs = [], {}
    sets = [("ncbi", ncbi), ("fid", fid)]
    recp = os.path.join(a.w, "fid", "records.tsv")
    if os.path.exists(recp):  # fid labels split by the species complex of the (first) FusariumID name
        cx = {}
        for r in read_tsv(recp):
            cx.setdefault(r["species"], r["complex"].split(" (")[0])
        groups = collections.defaultdict(dict)
        for g, l in fid.items():
            groups[cx.get(l.split(";")[0], "?")][g] = l
        for c, labs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            if len(labs) >= 20:
                sets.append((f"fid:{c}", labs))
    for name, labels in sets:
        rows = ani_gap(pairs, labels)
        if not rows:
            continue
        het = [h for *_x, h in rows]
        nn_ok = sum(c > h for _g, _l, c, h in rows)
        by_sp = collections.defaultdict(list)
        for _g, l, c, h in rows:
            by_sp[l].append((c, h))
        sp_gap = sum(min(c for c, _ in v) > max(h for _, h in v) for v in by_sp.values())
        gap_rows.append(
            dict(
                labels=name,
                genomes=len(rows),
                species=len(by_sp),
                nearest_is_same_species=f"{nn_ok / len(rows):.4f}",
                species_with_barcode_gap=f"{sp_gap / len(by_sp):.4f}",
                other_species_ani_q50=f"{q(het, 0.5):.4f}",
                other_species_ani_q90=f"{q(het, 0.9):.4f}",
                other_species_ani_q99=f"{q(het, 0.99):.4f}",
                other_species_ani_max=f"{max(het):.4f}",
                other_species_above_species_gate=f"{sum(h >= SPECIES_ANI for h in het) / len(het):.4f}",
                same_species_ani_q10=f"{q([c for *_x, c, _h in rows], 0.1):.4f}",
                same_species_ani_q50=f"{q([c for *_x, c, _h in rows], 0.5):.4f}",
            )
        )
    write_tsv(os.path.join(a.outdir, "ani_gap.tsv"), gap_rows)
    for x, y, ani in pairs:
        la, lb = fid.get(x), fid.get(y)
        if la and lb and la != lb:
            k = tuple(sorted((la, lb)))
            sp_pairs[k] = max(sp_pairs.get(k, 0.0), ani)
    top = sorted(sp_pairs.items(), key=lambda kv: -kv[1])[:60]
    write_tsv(
        os.path.join(a.outdir, "ani_species_pairs.tsv"),
        [dict(species_a=k[0], species_b=k[1], max_genome_ani=f"{v:.4f}") for k, v in top],
    )


if __name__ == "__main__":
    main()
