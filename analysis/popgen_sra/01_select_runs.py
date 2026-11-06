#!/usr/bin/env python3
"""01_select_runs.py : choose the SRA runs of a population set per species for popsketch.

From MATPredict's survey tables (wave_inputs/fungi_wgs_illumina.tsv and fungi_meta.tsv): paired Illumina WGS of the species, one run
per biosample (the smallest that still reaches --min-cov), mutant / evolved / cross / knockout looking records removed by their
titles, at most --cap per species taken round-robin across studies so no single study dominates. Writes runs.tsv per species with
the species' own genome size. Standard library only.   usage: 01_select_runs.py OUTDIR [--cap 300] [--min-cov 15]"""
import argparse, csv, os, re, collections

W = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/sra_reads_typing/wave_inputs/"
# species name as in the survey, slug, genome size in Mb
SPECIES = [("Zymoseptoria tritici", "z_tritici", 39), ("Cryptococcus neoformans", "c_neoformans", 19),
           ("Cercospora beticola", "c_beticola", 37), ("Parastagonospora nodorum", "p_nodorum", 37),
           ("Botrytis cinerea", "b_cinerea", 42), ("Neurospora discreta", "n_discreta", 37),
           ("Pyricularia oryzae", "p_oryzae", 42)]  # = Magnaporthe oryzae
BAD = re.compile(r"mutant|knock|delet|disrupt|evolv|\bko\b|crispr|rnai|tn-?seq|insertion|transform|overexpress|cross\b|progeny|offspring|"
                 r"\bF[12]\b|backcross|segregant|tetrad|lab[- ]?adapt|passage|ancestor|time ?point|generation|mutation accumulation|"
                 r"experimental evolution|\bT-DNA|\bATMT|heterokaryon|diploid|hybrid", re.I)

ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--cap", type=int, default=300); ap.add_argument("--min-cov", type=float, default=15)
a = ap.parse_args()
meta = {}
with open(W + "fungi_meta.tsv", newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        meta[r["run_accession"]] = r
by = collections.defaultdict(list)
with open(W + "fungi_wgs_illumina.tsv", newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        by[r["scientific_name"].split(" var.")[0].split(" f. sp.")[0]].append(r)
for name, slug, mb in SPECIES:
    rows = [r for k, v in by.items() if k == name for r in v]
    n_all = len({r["sample_accession"] for r in rows})
    ok = []
    for r in rows:
        if r["library_layout"] != "PAIRED" or float(r["base_count"] or 0) < a.min_cov * mb * 1e6: continue
        m = meta.get(r["run_accession"], {})
        text = " ".join(m.get(k, "") for k in ("study_title", "sample_title", "experiment_title"))
        if BAD.search(text): continue
        ok.append(r)
    best = {}
    for r in ok:  # one run per biosample: the smallest that passes
        s = r["sample_accession"]
        if s not in best or float(r["base_count"]) < float(best[s]["base_count"]): best[s] = r
    studies = collections.defaultdict(list)
    for r in sorted(best.values(), key=lambda r: r["run_accession"]): studies[r["study_accession"]].append(r)
    pick = []
    queues = [list(v) for _, v in sorted(studies.items(), key=lambda kv: -len(kv[1]))]
    while len(pick) < a.cap and any(queues):
        for q in queues:
            if q and len(pick) < a.cap: pick.append(q.pop(0))
    d = os.path.join(a.out, slug); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "runs.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["run_accession", "study_accession", "sample_accession", "species", "strain", "base_count", "read_count", "coverage", "center_name", "country", "genome_mb"])
        for r in pick:
            m = meta.get(r["run_accession"], {})
            w.writerow([r["run_accession"], r["study_accession"], r["sample_accession"], name, m.get("strain") or m.get("isolate") or m.get("sample_alias", ""),
                        r["base_count"], r["read_count"], round(float(r["base_count"]) / (mb * 1e6), 1), r["center_name"], m.get("country", ""), mb])
    print(f"{name}: {n_all} biosamples in the survey, {len(best)} pass (paired, >= {a.min_cov}x, no mutant-like title), {len(pick)} chosen from {len({r['study_accession'] for r in pick})} studies")
