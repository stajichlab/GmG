#!/usr/bin/env python3
"""32_bd_ploidy_from_ad.py : are the Bd samples that popsketch still flags `abundance_shoulder` varying ploidy or contamination?

At heterozygous SNP calls the fraction of reads carrying the alternate allele (from the VCF's AD field) shows the ploidy: a
diploid is centred on 0.5, a triploid has peaks near 1/3 and 2/3, a tetraploid near 1/4, 1/2 and 3/4. A contaminant gives a skew
to low fractions, not a regular series. Compares the flagged samples with unflagged controls of similar coverage.
Run where bcftools is on the PATH (HPCC: module load bcftools/1.19).   usage: 32_bd_ploidy_from_ad.py"""
import json, re, subprocess, sys
import numpy as np, pandas as pd

P = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/"
T = "/rhome/jstajich/projects/TE_fungi_copy/"
PD = "/bigdata/stajichlab/shared/projects/Population_Genomics/B_dendrobatidis/"
VCF = PD + "vcf/Bd_2.All.SNP.combined_selected.vcf.gz"
norm = lambda s: re.sub(r"^bd", "", re.sub(r"[^a-z0-9]", "", str(s).lower()))

d = json.load(open(P + "popsketch_runs/bd_full/popsketch.json"))
flagged = d["flags"].get("abundance_shoulder", [])
q = pd.read_csv(P + "popsketch_runs/bd_full/qc.tsv", sep="\t").set_index("id")
q["flags"] = q["flags"].fillna("")
ctrl_pool = q[(q.status == "included") & (q["flags"] == "") & q.kmer_coverage.between(20, 35)].index.tolist()
controls = ctrl_pool[:: max(1, len(ctrl_pool) // 12)][:12]
print(f"flagged: {flagged}\ncontrols ({len(controls)}): {controls}")

vs = [l.strip() for l in open(P + "popgen_bd/vcf_samples.txt")]
vmap = {}
for v in vs:
    vmap.setdefault(norm(v), []).append(v)
man = pd.read_csv(P + "popgen_bd/manifest.tsv", sep="\t").set_index("id")
meta = pd.read_csv(T + "results/bd/bd_metadata.tsv", sep="\t").set_index("strain")
tm = pd.read_csv(T + "results/bd/relatedness/bd_tree_map.tsv", sep="\t").set_index("strain")
sc = pd.read_csv(PD + "samples.csv"); srr2strain = dict(zip(sc.FileBase, sc.Strain))
def to_vcf(p):
    cands = [p, srr2strain.get(p)]
    if p in tm.index: cands.append(tm.tree_tip.get(p))
    if p in meta.index: cands += [meta.sheet_strain.get(p), meta.isolate.get(p)]
    for c in cands:
        if isinstance(c, str) and c and len(vmap.get(norm(c), [])) == 1:
            return vmap[norm(c)][0]
want = {p: to_vcf(p) for p in flagged + controls}
print("mapped to the VCF:", {p: v for p, v in want.items() if v}, "| not in the VCF:", [p for p, v in want.items() if not v])
names = [v for v in want.values() if v]
back = {v: p for p, v in want.items() if v}

cmd = f"bcftools view -s {','.join(names)} -m2 -M2 -v snps {VCF} | bcftools query -f '[%SAMPLE\\t%GT\\t%AD\\n]'"
fr = {n: [] for n in names}
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    s, gt, ad = line.rstrip("\n").split("\t")
    if gt not in ("0/1", "0|1", "1|0") or "," not in ad:
        continue
    try:
        r, a = (int(x) for x in ad.split(",")[:2])
    except ValueError:
        continue
    if r + a >= 15:
        fr[s].append(a / (r + a))
proc.wait()

rows = []
for n in names:
    f = np.array(fr[n])
    if len(f) < 200:
        rows.append(dict(id=back[n], flagged=back[n] in flagged, het_sites=len(f))); continue
    rows.append(dict(id=back[n], flagged=back[n] in flagged, het_sites=len(f), mean=f.mean(),
                     near_half=((f >= 0.42) & (f <= 0.58)).mean(), near_third=(((f >= 0.27) & (f < 0.40)) | ((f > 0.60) & (f <= 0.73))).mean(),
                     near_quarter=(((f >= 0.18) & (f < 0.27)) | ((f > 0.73) & (f <= 0.82))).mean(), low_tail=(f < 0.18).mean() + (f > 0.82).mean()))
r = pd.DataFrame(rows).set_index("id")
pd.set_option("display.width", 200)
print("\nshare of heterozygous calls (depth >= 15) by alt-allele fraction:")
print(r.round(3).to_string())
if "near_half" in r:
    print("\nmedians:"); print(r.groupby("flagged")[["het_sites", "near_half", "near_third", "near_quarter", "low_tail"]].median().round(3).to_string())
