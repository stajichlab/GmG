#!/usr/bin/env python3
"""29_check_popsketch.py : compare popsketch's flags on the real read sets with the problems found earlier by other means."""
import sys
import pandas as pd
R = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsketch_runs/"
suffix = sys.argv[1] if len(sys.argv) > 1 else "_full"
def load(ds):
    q = pd.read_csv(R + ds + suffix + "/qc.tsv", sep="\t").rename(columns={"flags": "flg"}).set_index("id"); q["flg"] = q["flg"].fillna(""); return q
q = load("cocci")
def show(title, ids):
    print(f"\n{title}")
    for i in ids:
        if i in q.index: print(f"  {i:12s} {q.status[i]:9s} {q.flg[i]}")
        else: print(f"  {i:12s} (not in run)")
show("not Coccidioides (Histoplasma, Lecanicillium x2, Fusarium)", ["NM_7898", "CA10", "CA12", "TX13"])
show("index-confirmed heavy mixtures (31-46%)", ["TX7", "M192", "NM_0071", "M158", "M216"])
show("B12526 (20% contaminant) and M223 (7%)", ["B12526", "M223"])
g = pd.read_csv("/rhome/jstajich/projects/TE_fungi_copy/results/cocci_gmg/out/results.tsv", sep="\t").set_index("query_id")
mix = g[g.qc_flags.fillna("").str.contains("mixture_suspected")].index
hit = [i for i in mix if i in q.index and q.flg[i] != ""]
print(f"\nindex-confirmed mixtures in the earlier gmg run: {len(mix)}; flagged by popsketch (any flag): {len(hit)}; with mixture_suspected: {sum('mixture_suspected' in q.flg[i] for i in mix if i in q.index)}")
print("  not flagged:", [i for i in mix if i in q.index and q.flg[i] == ""])
print("\nflag counts:", q.flg.str.split(";").explode().value_counts().drop("", errors="ignore").to_dict())
lowcov = q[q.flg.str.contains("low_coverage")]
print("low_coverage samples:", len(lowcov), "(python-era count was 14 at < 4x)")
other = q[(q.flg != "") & ~q.flg.str.contains("low_coverage") & ~q.index.isin(["NM_7898","CA10","CA12","TX13","TX7","M192","NM_0071","M158","M216","B12526","M223"]) & ~q.index.isin(mix)]
print("\nflagged for other reasons (not low coverage, not in the lists above):"); print(other[["status", "flg"]].to_string())
b = load("bd")
print("\n===== Bd")
half = pd.read_csv(R.replace("popsketch_runs/", "popgen_bd/") + "out_s200_qc/qc.tsv", sep="\t")
half_ids = set(half[half.status == "excluded"].id)
print("flag counts:", b.flg.str.split(";").explode().value_counts().drop("", errors="ignore").to_dict())
sl = set(b[b.flg.str.contains("size_low")].index)
rust = pd.read_csv("/bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_bd/rust_cutoff_bd_s200.tsv", sep="\t").set_index("id")
known = set(rust[rust.est_mb < 17].index)
print(f"half-content samples found earlier (est < 17 Mb with the Rust cutoff): {len(known)}; popsketch size_low: {len(sl)}; overlap {len(known & sl)}; missed {sorted(known - sl)}; extra {sorted(sl - known)[:10]}")
gen = pd.read_csv("/rhome/jstajich/projects/TE_fungi_copy/resources/panels/bd_species_audit.tsv", sep="\t").set_index("strain")
print("Sample_55930 (index-confirmed 6.3% mixture):", b.loc["Sample_55930", ["status", "flg"]].tolist() if "Sample_55930" in b.index else "n/a")
print("shoulder-flagged Bd samples:", sorted(b[b.flg.str.contains("abundance_shoulder")].index))
