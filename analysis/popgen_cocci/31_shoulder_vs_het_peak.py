#!/usr/bin/env python3
"""31_shoulder_vs_het_peak.py : can the shoulder statistic be made to ignore heterozygosity?

For each sample with k-mer coverage c >= 12 (from popsketch's qc.tsv and sketches): s = hashes at 0.1c-0.4c, h = hashes at
0.4c-0.65c (where a heterozygous peak sits), both as a fraction of the genome mass (hashes at >= 0.6c). A contaminant adds to
s without adding to h; heterozygosity adds to h and, through its lower tail, a roughly proportional amount to s. Reports the
robust z (cohort, coverage bins as in popsketch) of s, and of the residual of s after regressing it on h within the cohort.
Known: Cocci B12526 is a 20% contaminant at about 0.22c; the 11 Bd samples flagged `abundance_shoulder` were suspected to be
heterozygous.   usage: 31_shoulder_vs_het_peak.py"""
import os, struct
import numpy as np, pandas as pd

R = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/popsketch_runs/"


def abunds(path):
    with open(path, "rb") as f:
        f.read(8); k, scaled, seed, nl = struct.unpack("<IIQI", f.read(20)); f.read(nl)
        (n,) = struct.unpack("<Q", f.read(8)); f.seek(8 * n, 1); f.read(1)
        return np.frombuffer(f.read(4 * n), dtype="<u4")


def feats(ds):
    q = pd.read_csv(R + ds + "_full/qc.tsv", sep="\t").set_index("id")
    rows = []
    for sid, r in q.iterrows():
        c = r.kmer_coverage
        p = R + ds + "_full/sketches/" + sid + ".gsk"
        if not os.path.exists(p) or not np.isfinite(c) or c < 12 or r.status == "error":
            continue
        a = abunds(p)  # raw abundances, as popsketch's QC uses them: a contaminant sits below the abundance cutoff
        mass = (a >= 0.6 * c).sum()
        if mass < 1000:
            continue
        lo = max(3, int(0.1 * c))
        s = ((a >= lo) & (a < 0.4 * c)).sum() / mass
        h = ((a >= 0.4 * c) & (a < 0.65 * c)).sum() / mass
        rows.append(dict(id=sid, c=c, s=s, h=h, flags=str(r["flags"]) if isinstance(r["flags"], str) else ""))
    d = pd.DataFrame(rows).set_index("id")
    bins = pd.cut(d.c, [0, 15, 25, 40, 1000])
    # residual of s after a robust line s ~ a + b h within the cohort (Theil-Sen-like: median slope from random pairs)
    rng = np.random.default_rng(0)
    i, j = rng.integers(0, len(d), 4000), rng.integers(0, len(d), 4000)
    ok = np.abs(d.h.values[i] - d.h.values[j]) > 1e-3
    slope = np.median((d.s.values[i][ok] - d.s.values[j][ok]) / (d.h.values[i][ok] - d.h.values[j][ok]))
    icpt = np.median(d.s.values - slope * d.h.values)
    d["resid"] = d.s - (icpt + slope * d.h)
    for col in ("s", "resid"):
        z = pd.Series(np.nan, index=d.index)
        for b, idx in d.groupby(bins, observed=True).groups.items():
            v = d.loc[idx, col]
            if len(v) < 8: continue
            m = v.median(); sd = max(1.4826 * (v - m).abs().median(), 0.001)
            z.loc[idx] = (v - m) / sd
        d[col + "_z"] = z
    return d, slope, icpt


for ds in ("bd", "cocci"):
    d, slope, icpt = feats(ds)
    shoulder = d["flags"].str.contains("abundance_shoulder")
    print(f"\n===== {ds}: {len(d)} samples with c >= 12; slope of s on h = {slope:.3f}; flagged by the current rule: {int(shoulder.sum())}")
    print(f"  samples with s_z > 5: {int((d.s_z > 5).sum())}; with residual z > 5: {int((d.resid_z > 5).sum())}")
    print("  currently flagged samples, s_z vs residual z:")
    print(d[shoulder][["c", "s", "h", "s_z", "resid_z"]].round(3).to_string())
    if ds == "cocci":
        print("  known contaminants and other high-s samples:")
        print(d.loc[[x for x in ("B12526", "M223", "B13956") if x in d.index], ["c", "s", "h", "s_z", "resid_z"]].round(3).to_string())
        print("  highest residual z:"); print(d.sort_values("resid_z", ascending=False).head(6)[["c", "s", "h", "s_z", "resid_z", "flags"]].round(3).to_string())
