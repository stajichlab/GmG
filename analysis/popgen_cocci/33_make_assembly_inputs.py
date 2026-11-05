#!/usr/bin/env python3
"""33_make_assembly_inputs.py : sample sheets and labels for `popsketch run` on real BFD assemblies.

asm_cocci: the 22 Coccidioides reference assemblies plus 3 Histoplasma and 2 Aspergillus fumigatus assemblies as planted outgroups.
asm_asp:   every Aspergillus assembly in the BFD set (about 1,600), labelled by species: a real set above 560 samples."""
import os
import pandas as pd
B = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/"
O = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsketch_runs/"
m = pd.read_csv(B + "samples.csv")
m["path"] = B + "input_clean_genomes/" + m.ASMID + ".fa.gz"
m = m[[os.path.exists(p) for p in m.path]]
def write(name, df):
    with open(O + f"{name}_samples.tsv", "w") as f:
        f.write("id\tpath\n"); [f.write(f"{r.ASMID}\t{r.path}\n") for r in df.itertuples()]
    with open(O + f"{name}_labels.tsv", "w") as f:
        f.write("id\tspecies\tgenus\n"); [f.write(f"{r.ASMID}\t{r.SPECIES}\t{r.GENUS}\n") for r in df.itertuples()]
    print(name, len(df), "assemblies;", df.SPECIES.nunique(), "species")
cocci = m[m.GENUS == "Coccidioides"]
out = pd.concat([m[m.GENUS == "Histoplasma"].head(3), m[m.SPECIES == "Aspergillus fumigatus"].head(2)])
write("asm_cocci", pd.concat([cocci, out]))
write("asm_asp", m[m.GENUS == "Aspergillus"])
print("planted outgroups:", out.ASMID.tolist())
