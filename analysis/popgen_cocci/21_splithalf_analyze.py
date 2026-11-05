#!/usr/bin/env python3
"""21_splithalf_analyze.py : does a split-half noise floor predict where sketch distances stop carrying structure?

Floor = 1 - ANI between the two halves (alternate reads) of one sample, using the same abundance cutoff rule (gmg-core's,
via cutoff_probe) and the same containment ANI as the cohort matrices. Prediction under test (from the plan):
  (1) Bd lineages 2/3/4 should be mostly *below* the floor (they cannot be told apart), Bd between-group pairs above it;
  (2) among pairs below the floor, sketch distance should agree with SNP-PC distance much less than among pairs above it.
usage: 21_splithalf_analyze.py   (reads splithalf/*/sk/*.gsk, splithalf_cutoffs.tsv, the *_rc cohort ANI matrices)"""
import csv, os, struct, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr

P = "/bigdata/stajichlab/jstajich/projects/GmG_popgen/"


def read_gsk(path):
    with open(path, "rb") as f:
        assert f.read(8) == b"GMGSKT01"
        k, scaled, seed, nl = struct.unpack("<IIQI", f.read(20)); f.read(nl)
        (n,) = struct.unpack("<Q", f.read(8)); h = np.frombuffer(f.read(8 * n), dtype="<u8")
        has = f.read(1)[0]; a = np.frombuffer(f.read(4 * n), dtype="<u4") if has else None
    return k, h, a


cut = pd.read_csv(P + "splithalf_cutoffs.tsv", sep="\t").set_index("id").rust_min_abund
cfg = {"cocci": dict(ani=P + "popgen_cocci/final_s200_rc/ani_est.tsv", snp=P + "popgen_cocci/snp_pca_final_s200.tsv",
                     groups=P + "popgen_cocci/final_s200_rc/groups.tsv"),
       "bd": dict(ani=P + "popgen_bd/out_s200_qc_rc/ani_est.tsv", snp=P + "popgen_bd/snp_pca_vcf.tsv",
                  groups=None)}
bdlab = pd.read_csv(P + "popgen_bd/manifest.tsv", sep="\t").set_index("id").organism
for ds, c in cfg.items():
    d = f"{P}splithalf/{ds}/sk/"
    ids = sorted({f[:-6] for f in os.listdir(d) if f.endswith("_a.gsk")} & {f[:-6] for f in os.listdir(d) if f.endswith("_b.gsk")})
    rows = []
    for i in ids:
        k, ha, aa = read_gsk(d + i + "_a.gsk"); _, hb, ab = read_gsk(d + i + "_b.gsk")
        a_ = ha[aa >= cut.get(i + "_a", 2)]; b_ = hb[ab >= cut.get(i + "_b", 2)]
        inter = len(np.intersect1d(a_, b_, assume_unique=True)); cont = inter / max(1, min(len(a_), len(b_)))
        rows.append(dict(id=i, cut_a=cut.get(i + "_a", 2), cut_b=cut.get(i + "_b", 2), n_a=len(a_), n_b=len(b_), ani=cont ** (1 / k)))
    sh = pd.DataFrame(rows).set_index("id"); sh["dist"] = 1 - sh.ani
    print(f"\n===== {ds}: {len(sh)} samples, split-half distance (1-ANI)")
    print(sh[["cut_a", "n_a", "n_b", "dist"]].describe().loc[["min", "50%", "max"]].round(5).to_string())
    A = pd.read_csv(c["ani"], sep="\t"); A.index = A.columns
    D = 1 - A
    iu = np.triu_indices(len(D), 1); dv = D.values[iu]
    floor_med, floor_hi = sh.dist.median(), sh.dist.quantile(0.9)
    print(f"floor (median) {floor_med:.5f}, floor (90th pct) {floor_hi:.5f}; cohort pairs n={len(dv)}; "
          f"fraction of pairs below median floor {np.mean(dv < floor_med):.3f}, below 90th-pct floor {np.mean(dv < floor_hi):.3f}")
    if ds == "bd":
        coarse = {"Bd_L2": "L2-4", "Bd_L3": "L2-4", "Bd_L4": "L2-4", "Bd_L5": "L5/7", "Bd_L7": "L5/7", "Bd_L1": "L1", "Bd_L6": "L6", "Bd_L8": "L8"}
        g = bdlab.reindex(D.index).map(coarse)
        a = g.index[g == "L2-4"]; o = g.index[g.isin(["L1", "L6", "L8", "L5/7"])]
        w = D.loc[a, a].values[np.triu_indices(len(a), 1)]; b = D.loc[a, o].values.ravel()
        print(f"Bd within L2/3/4: median dist {np.median(w):.5f}, fraction below floor(median) {np.mean(w < floor_med):.3f} / (90th) {np.mean(w < floor_hi):.3f}")
        print(f"Bd L2/3/4 vs other groups: median dist {np.median(b):.5f}, fraction below floor(median) {np.mean(b < floor_med):.3f} / (90th) {np.mean(b < floor_hi):.3f}")
    else:
        g = pd.read_csv(c["groups"], sep="\t").set_index("id").organism_label.reindex(D.index)
        for sp in ("Coccidioides immitis", "Coccidioides posadasii"):
            a = g.index[g == sp]; w = D.loc[a, a].values[np.triu_indices(len(a), 1)]
            print(f"{sp} within-species: median dist {np.median(w):.5f}, fraction below floor(median) {np.mean(w < floor_med):.3f} / (90th) {np.mean(w < floor_hi):.3f}")
    # agreement with SNP-PC distance above vs below the floor
    S = pd.read_csv(c["snp"], sep="\t", index_col=0); sid = [i for i in D.index if i in S.index]
    Dm = D.loc[sid, sid].values; Sm = np.sqrt(((S.loc[sid].values[:, None, :5] - S.loc[sid].values[None, :, :5]) ** 2).sum(-1))
    iu2 = np.triu_indices(len(sid), 1); kd, sd = Dm[iu2], Sm[iu2]
    for name, fl in (("median floor", floor_med), ("90th-pct floor", floor_hi)):
        lo, hi = kd < fl, kd >= fl
        r_lo = spearmanr(kd[lo], sd[lo]).correlation if lo.sum() > 50 else float("nan")
        r_hi = spearmanr(kd[hi], sd[hi]).correlation if hi.sum() > 50 else float("nan")
        print(f"  Spearman(sketch dist, SNP-PC dist) {name}: below floor rho={r_lo:.3f} (n={lo.sum()}), above floor rho={r_hi:.3f} (n={hi.sum()})")
