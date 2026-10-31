#!/usr/bin/env python3
"""Score the long-read classifications of #44 (simulated sweep and real runs).

usage:
  long_reads_score.py collect-sim SIM_ROOT GENOMES.tsv OUT_QUERIES.tsv
      SIM_ROOT/<genome id>/p<k>/{measured.tsv, <tag>/results.tsv, <tag>/results.jsonl}; tag = <hidden level>_a<min-abund> (self_a2, species_a2, ...).
      One row per classified read set: genome, truth lineage, nominal error, measured error, depth, tag, call at every rank, status, ANI, containment,
      shared hashes, n_hashes, coverage_est, flags, the abundance histogram (counts of hashes seen 1..12 times, 13+, total bases, reads).
  long_reads_score.py collect-real RUNS.tsv REAL_ROOT OUT_QUERIES.tsv
      REAL_ROOT/<tag>/results.{tsv,jsonl}, the ids are run accessions; truth from RUNS.tsv (data/validation/long_reads_runs.tsv: lineage columns,
      measured error); same columns as above (depth = the coverage of the first 200 Mb of reads over the genome size).
  long_reads_score.py summarize QUERIES.tsv OUT_ACCURACY.tsv
      precision of the calls made, called share and recall over the attainable inputs per (tag, nominal error, depth, rank), medians of ANI and coverage_est.
Truth = the lineage of the genome the reads were taken from. Attainable: self hidden: genus if the genus has another genome in the reference, species if
the species has one; species hidden: genus if another species of the genus remains, species never; genus hidden: neither.
"""
import csv, json, os, sys, glob, statistics, collections

RANKS = ["phylum", "class", "order", "family", "genus", "species"]


def read_tsv(p):
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def hist_fields(rec):
    q = (rec or {}).get("qc", {})
    h = q.get("abundance_hist") or {}
    c = list(h.get("counts", [0] * 12)) + [h.get("over12", 0)]
    return {f"h{i+1}": c[i] for i in range(12)} | {"h13p": c[12], "median_kept": h.get("median_kept", ""), "total_bases": q.get("total_bases", ""),
                                                    "n_records": q.get("n_records", ""), "n_hashes": q.get("n_hashes", "")}


def load_dir(d):
    rows = {r["query_id"]: r for r in read_tsv(f"{d}/results.tsv")}
    jl = {}
    if os.path.exists(f"{d}/results.jsonl"):
        for line in open(f"{d}/results.jsonl"):
            try:
                r = json.loads(line)
            except ValueError:
                continue
            jl[r.get("query_id") or r.get("id")] = r
    return rows, jl


OUT_COLS = (["source", "genome", "species", "genus", "phylum", "tag", "hidden", "min_abund", "nominal_error", "measured_error", "depth", "status", "best_rank"]
            + [f"call_{r}" for r in RANKS] + ["tentative_genus", "ani", "containment_query", "shared", "coverage_est", "qc_flags", "withheld", "top_hit_species",
            "species_attainable", "genus_attainable"] + [f"truth_{r}" for r in RANKS]
            + [f"h{i}" for i in range(1, 13)] + ["h13p", "median_kept", "total_bases", "n_records", "n_hashes"])


def attain(g, hidden):
    so, go = int(g["species_others_in_ref"]), int(g["genus_others_in_ref"])
    if hidden == "self": return so > 0, go > 0
    if hidden == "species": return False, go - so > 0
    return False, False


def collect_sim(root, genomes_p, out):
    G = {g["id"]: g for g in read_tsv(genomes_p)}
    rows = []
    for gdir in sorted(glob.glob(f"{root}/*/")):
        gid = os.path.basename(gdir.rstrip("/"))
        if gid not in G: continue
        g = G[gid]
        for pdir in sorted(glob.glob(f"{gdir}p*/")):
            meas = {m["label"]: m for m in read_tsv(f"{pdir}measured.tsv")} if os.path.exists(f"{pdir}measured.tsv") else {}
            for tdir in sorted(glob.glob(f"{pdir}*_a*/")):
                tag = os.path.basename(tdir.rstrip("/"))
                if not os.path.exists(f"{tdir}results.tsv"): continue
                hidden, ma = tag.split("_a")
                res, jl = load_dir(tdir.rstrip("/"))
                for qid, r in res.items():
                    e, d = qid[1:].split("_d")
                    sa, ga = attain(g, hidden)
                    row = {"source": "sim", "genome": gid, "species": g["species"], "genus": g["genus"], "phylum": g["phylum"], "tag": tag, "hidden": hidden,
                           "min_abund": ma, "nominal_error": float(e) * 100, "measured_error": (meas.get(f"e{e}") or {}).get("error_pct", ""), "depth": d,
                           "species_attainable": int(sa), "genus_attainable": int(ga)}
                    for k in ["status", "best_rank", "tentative_genus", "ani", "containment_query", "shared", "coverage_est", "qc_flags", "withheld", "top_hit_species"] + [f"call_{x}" for x in RANKS]:
                        row[k] = r.get(k, "")
                    for x in RANKS: row[f"truth_{x}"] = g[x]
                    row.update(hist_fields(jl.get(qid)))
                    rows.append(row)
    write(out, rows)
    print(len(rows), "rows", file=sys.stderr)


def collect_real(runs_p, root, out):
    R = {r["run"]: r for r in read_tsv(runs_p)}
    rows = []
    for tdir in sorted(glob.glob(f"{root}/*_a*/")):
        tag = os.path.basename(tdir.rstrip("/"))
        if not os.path.exists(f"{tdir}results.tsv"): continue
        hidden, ma = tag.split("_a")
        res, jl = load_dir(tdir.rstrip("/"))
        for qid, r in res.items():
            if qid not in R: continue
            g = R[qid]
            sa, ga = attain(g, hidden)
            cov = float(g["base_count_used"]) / float(g["genome_size"]) if g.get("base_count_used") else ""
            row = {"source": "real", "genome": g["bfd_id"], "species": g["species"], "genus": g["genus"], "phylum": g["phylum"], "tag": tag, "hidden": hidden,
                   "min_abund": ma, "nominal_error": g.get("measured_error_pct", ""), "measured_error": g.get("measured_error_pct", ""), "depth": f"{cov:.3g}" if cov != "" else "",
                   "species_attainable": int(sa), "genus_attainable": int(ga)}
            for k in ["status", "best_rank", "tentative_genus", "ani", "containment_query", "shared", "coverage_est", "qc_flags", "withheld", "top_hit_species"] + [f"call_{x}" for x in RANKS]:
                row[k] = r.get(k, "")
            for x in RANKS: row[f"truth_{x}"] = g[x]
            row["run"] = qid
            row.update(hist_fields(jl.get(qid)))
            rows.append(row)
    write(out, rows, extra=["run"])
    print(len(rows), "rows", file=sys.stderr)


def write(out, rows, extra=()):
    with open(out, "w") as f:
        w = csv.DictWriter(f, OUT_COLS + list(extra), delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)


def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def med(v):
    v = [x for x in v if x is not None]
    return statistics.median(v) if v else None


def summarize(qp, out):
    Q = read_tsv(qp)
    by = collections.defaultdict(list)
    for r in Q: by[(r["tag"], r["nominal_error"], r["depth"])].append(r)
    cols = ["tag", "nominal_error_pct", "depth", "n", "rank", "called", "precision", "n_called", "attainable", "recall", "ok_share", "no_match_share", "median_ani",
            "median_coverage_est", "median_cov_est_over_depth", "median_n_hashes"]
    rows = []
    for (tag, e, d), rs in sorted(by.items(), key=lambda kv: (kv[0][0], float(kv[0][1]), float(kv[0][2]))):
        n = len(rs)
        base = {"tag": tag, "nominal_error_pct": e, "depth": d, "n": n,
                "ok_share": sum(r["status"] == "ok" for r in rs) / n, "no_match_share": sum(r["status"] == "no_match" for r in rs) / n,
                "median_ani": med([fnum(r["ani"]) for r in rs]), "median_coverage_est": med([fnum(r["coverage_est"]) for r in rs]),
                "median_cov_est_over_depth": med([fnum(r["coverage_est"]) / float(d) if fnum(r["coverage_est"]) is not None else None for r in rs]),
                "median_n_hashes": med([fnum(r["n_hashes"]) for r in rs])}
        for rk in ("phylum", "family", "genus", "species"):
            called = [r for r in rs if r[f"call_{rk}"]]
            right = [r for r in called if r[f"call_{rk}"] == r[f"truth_{rk}"]]
            att = [r for r in rs if rk not in ("genus", "species") or int(r[f"{rk}_attainable"])]
            att_right = [r for r in att if r[f"call_{rk}"] and r[f"call_{rk}"] == r[f"truth_{rk}"]]
            rows.append(base | {"rank": rk, "called": len(called) / n, "precision": len(right) / len(called) if called else "", "n_called": len(called),
                                "attainable": len(att), "recall": len(att_right) / len(att) if att else ""})
    with open(out, "w") as f:
        w = csv.DictWriter(f, cols, delimiter="\t", lineterminator="\n"); w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(len(rows), "rows", file=sys.stderr)


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "collect-sim": collect_sim(*sys.argv[2:5])
    elif c == "collect-real": collect_real(*sys.argv[2:5])
    elif c == "summarize": summarize(*sys.argv[2:4])
    else: sys.exit(__doc__)
