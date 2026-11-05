#!/usr/bin/env python3
"""05_tree_check.py CLEAN_DIR [drop ids...] : rebuild NJ from mash_dist.tsv without the dropped ids, root on the species
split and count foreign tips inside each species' MRCA clade (0 = monophyletic). Run with and without label fixes."""
import sys, os
import numpy as np, pandas as pd
from io import StringIO
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from importlib import import_module
nj = import_module("02_popdist").nj
from Bio import Phylo
d, drop = sys.argv[1], set(sys.argv[2:])
D = pd.read_csv(f"{d}/mash_dist.tsv", sep="\t"); D.index = D.columns
keep = [i for i in D.columns if i not in drop]
D = D.loc[keep, keep]
lc = pd.read_csv(f"{d}/label_check.tsv", sep="\t").set_index("id")
for relabel in (False, True):
    L = {i: (lc.closest_label[i] if relabel and lc.conflict[i] else lc.label[i]) for i in keep}
    t = Phylo.read(StringIO(nj(D.values, keep)), "newick")
    for out in ("Coccidioides immitis", "Coccidioides posadasii"):
        pass
    t.root_with_outgroup(*[c for c in t.get_terminals() if L[c.name] == "Coccidioides immitis"])
    res = []
    for l in ("Coccidioides immitis", "Coccidioides posadasii"):
        tips = [c for c in t.get_terminals() if L[c.name] == l]
        foreign = [c.name for c in t.common_ancestor(tips).get_terminals() if L[c.name] != l]
        res.append(f"{l.split()[1]}: n={len(tips)} foreign_in_MRCA={len(foreign)} {foreign[:6]}")
    print(d, "labels corrected" if relabel else "labels as given", "|", " | ".join(res))
