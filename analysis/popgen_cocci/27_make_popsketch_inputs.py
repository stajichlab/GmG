#!/usr/bin/env python3
"""27_make_popsketch_inputs.py : sample sheets and label tables for `popsketch run` from the two manifests (Cocci, Bd).

Cocci labels: species (from samples.csv Organism), center. Bd labels: lineage (blank when unassigned), gmg_verdict."""
import csv, os, sys
import pandas as pd
P = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/"
OUT = P + "popsketch_runs/"
os.makedirs(OUT, exist_ok=True)

def sheet(manifest, path):
    m = pd.read_csv(manifest, sep="\t")
    n = 0
    with open(path, "w") as f:
        f.write("id\tpath\n")
        for _, r in m.iterrows():
            for col in ("r1", "r2"):
                p = r[col]
                if isinstance(p, str) and os.path.exists(p):
                    f.write(f"{r['id']}\t{p}\n"); n += 1
    return m, n

m, n = sheet(P + "popgen_cocci/manifest.tsv", OUT + "cocci_samples.tsv")
sc = pd.read_csv("/bigdata/stajichlab/shared/projects/Onygenales/Popgen/2025_All_Cocci/Assembly/samples.csv").rename(columns={"Strain": "id"}).set_index("id")
with open(OUT + "cocci_labels.tsv", "w") as f:
    f.write("id\tspecies\tcenter\n")
    for _, r in m.iterrows():
        c = sc.Center.get(r["id"]) if r["id"] in sc.index else None
        f.write(f"{r['id']}\t{r['organism']}\t{c if isinstance(c, str) else ''}\n")
print("cocci", len(m), "samples,", n, "files")
m, n = sheet(P + "popgen_bd/manifest.tsv", OUT + "bd_samples.tsv")
with open(OUT + "bd_labels.tsv", "w") as f:
    f.write("id\tlineage\tgmg_verdict\n")
    for _, r in m.iterrows():
        lin = "" if r["organism"] == "Bd_unassigned" else r["organism"]
        f.write(f"{r['id']}\t{lin}\t{r['notes'] if isinstance(r['notes'], str) else ''}\n")
print("bd", len(m), "samples,", n, "files")
