#!/usr/bin/env python3
"""Prepare the external validation sets (#40): a directory of symlinks per dataset and a labels table, on the HPCC.

usage: prep_external.py OUT_ROOT TAXDUMP_DIR       (OUT_ROOT/<set>/in/*.fa[.gz] and OUT_ROOT/<set>/labels.tsv)
labels.tsv columns: sample, name (provider name, may be empty), note (e.g. sample description such as metagenome).
Sets: kasson (KassonLab project genomes/ + samples.csv Species), antarctica (Euk MAGs, no names), ccfee_aaftf and
ccfee_shovill (the same strains assembled two ways; samples.csv Phylum), zygolife (aaftf/sort; NCBI taxid of the sample sheets).
Read-only on the source directories.
"""
import csv, glob, os, sys

P = "/bigdata/stajichlab/shared/projects"
out, td = sys.argv[1:3]


def link_dir(name):
    d = f"{out}/{name}/in"
    os.makedirs(d, exist_ok=True)
    for f in glob.glob(d + "/*"):
        os.unlink(f)
    return d


def write_labels(name, rows):
    with open(f"{out}/{name}/labels.tsv", "w") as o:
        w = csv.writer(o, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "name", "note"])
        w.writerows(rows)
    print(name, len(rows))


def pick(base, stem):
    """Prefer the plain FASTA over its .gz copy."""
    for ext in (".fasta", ".fasta.gz"):
        if os.path.exists(f"{base}/{stem}{ext}"):
            return f"{base}/{stem}{ext}", ext
    return None, None


# ---- KassonLab: per project, genomes/<strain>.AAFTF.fasta; names from the project's samples.csv.
# Some genomes/ directories also hold copies of another project's assemblies (Apophysomyces/ has Verticillium strains): a stem
# found in several projects is kept only where the project's own samples.csv lists it.
d = link_dir("kasson")
cand = {}   # (project, stem) -> (source, ext, species, note)
for proj in ["Apophysomyces", "Cryphonectria", "Entomophthorales", "Nectriaceae", "OphiocordycepsSymbiont", "Verticillium"]:
    pd = f"{P}/KassonLab/{proj}"
    sheet = {}
    for r in csv.DictReader(open(f"{pd}/samples.csv")):
        note = "; ".join(x for x in (r.get("Description", ""), r.get("AssemblyFocus", "")) if x)
        for k in ("Strain", "SeqStrain", "InternalID", "Sample", "Locus"):
            if r.get(k):
                sheet.setdefault(r[k].strip(), (r["Species"].strip(), note))
    for f in sorted(glob.glob(f"{pd}/genomes/*.AAFTF.fasta*")):
        stem = os.path.basename(f).split(".AAFTF.fasta")[0]
        if (proj, stem) in cand:
            continue
        src, ext = pick(f"{pd}/genomes", stem + ".AAFTF")
        sp, note = sheet.get(stem, ("", "not in this project's samples.csv"))
        cand[(proj, stem)] = (src, ext, sp, note)
labelled = {stem for (proj, stem), v in cand.items() if v[2]}
rows = []
for (proj, stem), (src, ext, sp, note) in cand.items():
    if not sp and stem in labelled:
        continue
    sid = f"{proj}__{stem}"
    os.symlink(src, f"{d}/{sid}.fa" + (".gz" if ext.endswith(".gz") else ""))
    rows.append([sid, sp, note])
write_labels("kasson", rows)

# ---- Antarctica eukaryotic MAGs: no names, classification only
d = link_dir("antarctica")
rows = []
for f in sorted(glob.glob(f"{P}/Antarctica/Euk_MAGs_Claudio/fasta/*")):
    b = os.path.basename(f)
    sid = b.split(".fa")[0]
    os.symlink(os.path.realpath(f), f"{d}/{b}")
    rows.append([sid, "", ""])
write_labels("antarctica", rows)

# ---- CCFEE: same strains, AAFTF vs shovill assemblies; samples.csv gives Phylum (species names are generic)
cc = f"{P}/Antarctica/Tuscia_LM_CCFEE_202205/assembly"
phy = {r["Strain"]: r["Phylum"] for r in csv.DictReader(open(f"{cc}/samples.csv"))}
for tag, ext in (("aaftf", "AAFTF"), ("shovill", "shovill")):
    d = link_dir(f"ccfee_{tag}")
    rows = []
    for f in sorted(glob.glob(f"{cc}/genomes/*.{ext}.fasta")):
        stem = os.path.basename(f).split(".")[0]
        os.symlink(f, f"{d}/{stem}.fa")
        key = stem.replace("MNA-", "")
        ph = next((v for k, v in phy.items() if k.replace("MNA-", "") == key or k == stem), "")
        rows.append([stem, "", f"phylum={ph}" if ph else ""])
    write_labels(f"ccfee_{tag}", rows)

# ---- ZyGoLife aaftf/sort: strain -> NCBI taxid in the sample sheets, name = scientific name of the taxid
za = f"{P}/ZyGoLife/LCG/Annotation_2026"
tax = {}
for f in sorted(glob.glob(f"{za}/samples*.csv")):
    for r in csv.DictReader(open(f)):
        if r.get("taxid"):
            tax.setdefault(r["sample"], r["taxid"])
have = sorted(x.replace(".sorted.fasta.gz", "") for x in os.listdir(f"{za}/results/aaftf/sort"))
want = {tax[s] for s in have if s in tax}
sci = {}
for l in open(f"{td}/names.dmp"):
    p = [x.strip() for x in l.split("|")]
    if p[3] == "scientific name" and p[0] in want:
        sci[p[0]] = p[1]
d = link_dir("zygolife")
rows = []
for s in have:
    os.symlink(f"{za}/results/aaftf/sort/{s}.sorted.fasta.gz", f"{d}/{s}.fa.gz")
    rows.append([s, sci.get(tax.get(s, ""), ""), f"taxid={tax.get(s, '')}"])
write_labels("zygolife", rows)
