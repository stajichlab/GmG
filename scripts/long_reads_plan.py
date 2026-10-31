#!/usr/bin/env python3
"""Pick public ONT / PacBio runs of fungi whose genome is in the reference, for the long-reads validation (#44).

usage:
  long_reads_plan.py candidates ENA_ONT.tsv ENA_PACBIO.tsv ASSEMBLY_SUMMARY.txt META_FULL.tsv META_FULL.tsv EXCLUSIONS.txt OUT.tsv [PER_GUESS=...] [SEED=44]
  long_reads_plan.py select CANDIDATES.tsv PILOT.tsv OUT.tsv [N=20]

ENA_*.tsv are `read_run` portal downloads (tax_tree(4751) AND instrument_platform=OXFORD_NANOPORE|PACBIO_SMRT AND library_source=GENOMIC; fields in
the header of the files: run_accession sample_accession secondary_sample_accession study_accession scientific_name instrument_platform instrument_model
library_strategy library_source library_layout read_count base_count fastq_ftp fastq_bytes center_name first_public experiment_title).

`candidates`: a run is kept when its BioSample is the BioSample of a genome of the reference (the full index, so that other genomes of the species stay when this one is hidden), the assembly is chromosome level or complete (NCBI assembly_level), the species is named the same way by the ENA
record, the NCBI organism name and the reference label (a binomial, not "sp."), the genome is not on the exclusion list, the run is single-ended with at
least 250 Mb of bases (the first 200 Mb must exist) and the genome is 2-60 Mb. The chemistry is only a GUESS from the metadata (year, instrument model, title);
candidates are stratified by the guess, at most one run per species and guess, at most 3 per genus. Writes OUT.tsv: the list that is pilot-measured
(scripts/long_reads_measure.py on the first 4 MB of each file) before anything is downloaded in full.
`select`: from the pilot's measured identity (mappy against the run's own assembly) take N runs spread over the four measured classes
(<0.5% HiFi-like, 0.5-3% R10.4-like, 3-11% R9.4-like, >=11% CLR-like), at most 2 runs per species and genus, as many species as possible.
"""
import csv, random, re, sys, collections, os

binom = lambda s: " ".join(re.sub(r"[\[\]]", "", s).split()[:2])


def good_name(s):
    w = binom(s).split()
    return "/" not in s and len(w) == 2 and w[0][:1].isupper() and w[1][:1].islower() and w[1] not in ("sp.", "sp", "cf.", "aff.") and "." not in w[1]


def guess(r):
    """chemistry guess from the metadata only"""
    year = int(r["first_public"][:4]); m = r["instrument_model"]; t = (r["experiment_title"] or "").lower()
    ml = int(r["base_count"]) / max(int(r["read_count"]), 1)
    if r["instrument_platform"] == "PACBIO_SMRT":
        if m in ("Revio", "Sequel IIe") or "hifi" in t or "ccs" in t or "circular consensus" in t: return "HiFi"
        if m == "Sequel II" and ml >= 8000 and year >= 2022: return "HiFi"
        return "CLR"
    if "r10" in t or "lsk114" in t or "q20" in t: return "ONT_R10"
    if year <= 2021: return "ONT_R9_old"
    if year >= 2024: return "ONT_R10"
    return "ONT_mid"


def candidates(a):
    ont, pb, asm_p, meta_p, reps_p, excl_p, out = a[:7]
    per_guess = int(a[7]) if len(a) > 7 else 16
    rng = random.Random(int(a[8]) if len(a) > 8 else 44)
    excl = {l.strip() for l in open(excl_p) if l.strip() and not l.startswith("#")}
    meta = {r["id"]: r for r in csv.DictReader(open(meta_p), delimiter="\t")}
    reps = {r["id"]: r for r in csv.DictReader(open(reps_p), delimiter="\t")}
    lines = open(asm_p).read().splitlines()
    hi = next(i for i, l in enumerate(lines) if l.startswith("#assembly_accession"))
    head = lines[hi].lstrip("#").split("\t")
    asm = {}
    for row in csv.reader(lines[hi + 1:], delimiter="\t"):
        d = dict(zip(head, row)); asm[d["assembly_accession"]] = d
    n_sp = collections.Counter(m["species"] for m in reps.values()); n_ge = collections.Counter(m["genus"] for m in reps.values())
    by_bs = collections.defaultdict(list)
    for gid, m in reps.items():
        acc = re.match(r"^(GC[AF]_\d+\.\d+)", gid)
        if not acc or gid in excl or acc.group(1) not in asm: continue
        A = asm[acc.group(1)]
        if A["assembly_level"] not in ("Complete Genome", "Chromosome") or not A["biosample"]: continue
        if not (good_name(m["species"]) and good_name(A["organism_name"])) or binom(m["species"]) != binom(A["organism_name"]): continue
        size = int(A["genome_size"] or 0)
        if not 2e6 <= size <= 60e6: continue
        by_bs[A["biosample"]].append((gid, A, size))
    rows = []
    for path in (ont, pb):
        for r in csv.DictReader(open(path), delimiter="\t"):
            hit = by_bs.get(r["sample_accession"]) or by_bs.get(r["secondary_sample_accession"])
            if not hit or r["library_layout"] != "SINGLE" or not r["read_count"] or not r["base_count"] or int(r["base_count"]) < 250e6: continue
            if not good_name(r["scientific_name"]) or not r["fastq_ftp"]: continue
            gid, A, size = hit[0]
            m = reps[gid]
            if binom(r["scientific_name"]) != binom(m["species"]): continue
            if ";" in r["fastq_ftp"] and "_1.fastq" in r["fastq_ftp"]: pass     # the first file is the single-end one in ENA's naming
            url = r["fastq_ftp"].split(";")[0]
            rows.append({"run": r["run_accession"], "platform": r["instrument_platform"], "model": r["instrument_model"], "year": r["first_public"][:4],
                         "chem_guess": guess(r), "species": m["species"], "genus": m["genus"], "phylum": m["phylum"], "class": m["class"],
                         "biosample": A["biosample"], "bfd_id": gid, "genome_size": size, "replicons": A["replicon_count"],
                         "species_others_in_ref": n_sp[m["species"]] - 1, "genus_others_in_ref": n_ge[m["genus"]] - 1,
                         "read_count": r["read_count"], "base_count": r["base_count"], "mean_len": int(int(r["base_count"]) / int(r["read_count"])),
                         "fastq_bytes": r["fastq_bytes"].split(";")[0], "center": r["center_name"], "title": (r["experiment_title"] or "")[:80], "url": url})
    rng.shuffle(rows)
    chosen, seen_sg, per_genus = [], set(), collections.Counter()
    for g in ("HiFi", "CLR", "ONT_R10", "ONT_mid", "ONT_R9_old"):
        n = 0
        for r in sorted((x for x in rows if x["chem_guess"] == g), key=lambda x: (x["species"] in {c["species"] for c in chosen}, x["species_others_in_ref"] == 0, rng.random())):
            if (r["species"], g) in seen_sg or per_genus[r["genus"]] >= 3 or r["mean_len"] < 2000: continue
            seen_sg.add((r["species"], g)); per_genus[r["genus"]] += 1; chosen.append(r); n += 1
            if n >= per_guess: break
    cols = list(chosen[0])
    with open(out, "w") as f:
        w = csv.DictWriter(f, cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(chosen)
    print(len(rows), "eligible runs;", len(chosen), "candidates", dict(collections.Counter(c["chem_guess"] for c in chosen)),
          len({c["species"] for c in chosen}), "species", file=sys.stderr)


def klass(err):
    return "HiFi" if err < 0.5 else "R10" if err < 3 else "R9" if err < 11 else "CLR"


def select(a):
    cand_p, pilot_p, out = a[:3]
    n = int(a[3]) if len(a) > 3 else 20
    rng = random.Random(44)
    cand = {r["run"]: r for r in csv.DictReader(open(cand_p), delimiter="\t")}
    rows = []
    for p in csv.DictReader(open(pilot_p), delimiter="\t"):
        if p["run"] in cand and p["error_pct"] not in ("", "NA"):
            r = dict(cand[p["run"]]); r["pilot_error_pct"] = p["error_pct"]; r["pilot_reads"] = p["n_aligned"]; r["pilot_median_len"] = p["median_len"]
            r["measured_class"] = klass(float(p["error_pct"])); rows.append(r)
    rng.shuffle(rows)
    quota = {"HiFi": 5, "R10": 5, "R9": 5, "CLR": 5}
    chosen, sp, ge = [], collections.Counter(), collections.Counter()
    for rnd in range(2):
        for k in quota:
            for r in rows:
                if r in chosen or r["measured_class"] != k or sp[r["species"]] > rnd or ge[r["genus"]] > rnd or sum(c["measured_class"] == k for c in chosen) >= quota[k]: continue
                if int(r["pilot_reads"]) < 50: continue
                chosen.append(r); sp[r["species"]] += 1; ge[r["genus"]] += 1
    chosen = chosen[:n]
    with open(out, "w") as f:
        w = csv.DictWriter(f, list(chosen[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(chosen)
    print(len(chosen), "runs;", dict(collections.Counter(c["measured_class"] for c in chosen)), len({c["species"] for c in chosen}), "species;",
          sum(int(c["fastq_bytes"]) for c in chosen) / 1e9, "GB in the full files (only the first 200 Mb of each is read)", file=sys.stderr)


def simgenomes(a):
    """Part 1: 24 chromosome-level genomes of the full index, over 5 phyla and as many classes as possible, 2.5-55 Mb, with a congener in the reference"""
    asm_p, reps_p, excl_p, out = a[:4]
    rng = random.Random(int(a[4]) if len(a) > 4 else 44)
    excl = {l.strip() for l in open(excl_p) if l.strip() and not l.startswith("#")}
    reps = {r["id"]: r for r in csv.DictReader(open(reps_p), delimiter="\t")}
    lines = open(asm_p).read().splitlines()
    hi = next(i for i, l in enumerate(lines) if l.startswith("#assembly_accession"))
    head = lines[hi].lstrip("#").split("\t")
    asm = {}
    for row in csv.reader(lines[hi + 1:], delimiter="\t"):
        d = dict(zip(head, row)); asm[d["assembly_accession"]] = d
    n_sp = collections.Counter(m["species"] for m in reps.values()); n_ge = collections.Counter(m["genus"] for m in reps.values())
    pool = collections.defaultdict(list)
    for gid, m in reps.items():
        acc = re.match(r"^(GC[AF]_\d+\.\d+)", gid)
        if not acc or gid in excl or acc.group(1) not in asm: continue
        A = asm[acc.group(1)]
        if A["assembly_level"] not in ("Complete Genome", "Chromosome"): continue
        if not (good_name(m["species"]) and good_name(A["organism_name"])) or binom(m["species"]) != binom(A["organism_name"]): continue
        size = int(A["genome_size"] or 0)
        if not 2.5e6 <= size <= 55e6 or n_ge[m["genus"]] < 2: continue
        pool[(m["phylum"], m["class"])].append({"id": gid, "species": m["species"], "genus": m["genus"], "phylum": m["phylum"], "class": m["class"],
            "order": m["order"], "family": m["family"], "genome_size": size, "replicons": A["replicon_count"],
            "species_others_in_ref": n_sp[m["species"]] - 1, "genus_others_in_ref": n_ge[m["genus"]] - 1})
    quota = {"Microsporidia": 2, "Mucoromycota": 1, "Mortierellomycota": 1, "Basidiomycota": 8, "Ascomycota": 12}
    chosen = []
    for ph, q in quota.items():
        classes = [k for k in pool if k[0] == ph]
        rng.shuffle(classes)
        n = 0
        while n < q and any(pool[k] for k in classes):
            for k in classes:
                if n >= q or not pool[k]: continue
                # half of the genomes: a conspecific rep exists (species recall is then attainable)
                cands = sorted(pool[k], key=lambda r: (-(r["species_others_in_ref"] > 0) if n % 2 == 0 else 0, rng.random()))
                r = cands[0]; pool[k].remove(r); chosen.append(r); n += 1
    with open(out, "w") as f:
        w = csv.DictWriter(f, list(chosen[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(chosen)
    print(len(chosen), "genomes;", len({c["class"] for c in chosen}), "classes;", sum(c["species_others_in_ref"] > 0 for c in chosen), "with a conspecific rep;",
          sum(c["genome_size"] for c in chosen) / 1e6, "Mb", file=sys.stderr)


if __name__ == "__main__":
    {"candidates": candidates, "select": select, "simgenomes": simgenomes}[sys.argv[1]](sys.argv[2:])
