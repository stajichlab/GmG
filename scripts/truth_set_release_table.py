#!/usr/bin/env python3
"""Truth-set tables of one run variant (a re-measurement on a later release), same cells as truth_set_tables.py (#39).

usage: truth_set_release_table.py SCORE_DIR VARIANT OUT.md      (SCORE_DIR has MODE.VARIANT.raw.refseq.tsv from validate_run.py)
Cells are called fraction / precision [95% Wilson] (calls made)."""
import importlib.util, os, sys, types

sd, variant, out = sys.argv[1:4]
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "truth_set_tables.py")).read().split("md = []")[0]
mod = types.ModuleType("tst")
sys.argv = ["x", sd, "/dev/null", "/dev/null"]
exec(src, mod.__dict__)
T = mod.table
M3 = [("loo", "leave-one-out"), ("los", "leave-species-out"), ("log", "leave-genus-out"), ("indep", "outside the index (21)")]
md = ["**All inputs**\n\n" + T(M3, variant=variant, kind="raw")]
for grp, mode, title in [("species_in_bfd=True", "loo", "leave-one-out, species still in the reference (answerable)"),
                         ("species_in_bfd=False", "loo", "leave-one-out, species a singleton in the reference (not answerable)"),
                         ("evidence=type_material", "loo", "leave-one-out, type material only"),
                         ("evidence=reference_genome", "loo", "leave-one-out, reference genomes (not type)"),
                         ("genus_in_bfd=True", "los", "leave-species-out, genus still in the reference"),
                         ("genus_in_bfd=False", "los", "leave-species-out, genus is a singleton (also absent)")]:
    t = T([(mode, title)], group=grp, variant=variant, kind="raw")
    if t.count("\n") > 2:
        md.append(f"**{title}**\n\n" + t)
open(out, "w").write("\n".join(md))
