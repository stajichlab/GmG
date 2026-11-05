#!/usr/bin/env python3
"""22_abundance_spectrum.py : can the k-mer abundance spectrum alone (no index) flag mixed samples?

Per sample (scaled-200 .gsk abundances; main coverage c and cutoff from the gmg-core probe):
  shoulder   = hashes with abundance in [lo, hi) as a fraction of the genome mass, lo = max(3, 0.1c), hi = 0.6c
  excess     = (observed - Poisson(c) expectation) in [lo, hi), as a fraction of the genome mass
  size_z     = robust z of the Rust-cutoff estimated size (cohort-relative size test, for comparison)
Both spectrum statistics are standardised against the cohort (median/MAD within coverage bins), because real coverage is
overdispersed. Truth: `gmg run --index` mixture results (second_explained_fraction); no label = not flagged.
Samples with c < MINC are 'not assessable' for the spectrum (no room between errors and the peak).

usage: 22_abundance_spectrum.py DATASET SKETCH_DIR PROBE.tsv GMG_RESULTS.tsv"""
import os, struct, sys
import numpy as np, pandas as pd
from scipy.stats import poisson
from sklearn.metrics import roc_auc_score

MINC = 6.0


def read_abunds(path):
    with open(path, "rb") as f:
        assert f.read(8) == b"GMGSKT01"
        k, scaled, seed, nl = struct.unpack("<IIQI", f.read(20)); f.read(nl)
        (n,) = struct.unpack("<Q", f.read(8)); f.seek(8 * n, 1)
        has = f.read(1)[0]
        return np.frombuffer(f.read(4 * n), dtype="<u4") if has else None


def features(a, c):
    h = np.bincount(np.minimum(a, 4000), minlength=4001).astype(float)       # h[x] = hashes with abundance x
    hi = int(np.ceil(0.6 * c)); lo = max(3, int(np.floor(0.1 * c)))
    if hi <= lo + 1:
        return np.nan, np.nan
    M = h[hi:].sum()                                                          # genome mass at/above 0.6c
    tail = poisson.sf(hi - 1, c)                                              # P(a >= hi | Poisson(c))
    if M <= 0 or tail <= 0:
        return np.nan, np.nan
    G = M / tail
    obs = h[lo:hi].sum(); exp = G * (poisson.cdf(hi - 1, c) - poisson.cdf(lo - 1, c))
    return obs / M, (obs - exp) / M


def robust_z(x, c):
    z = pd.Series(np.nan, index=x.index)
    bins = pd.cut(c, [0, 6, 10, 15, 25, 40, 1000])
    for b, idx in x.groupby(bins, observed=True).groups.items():
        v = x.loc[idx].dropna()
        if len(v) < 8: continue
        m = v.median(); s = 1.4826 * (v - m).abs().median() or 1e-6
        z.loc[idx] = (x.loc[idx] - m) / s
    return z


ds, skdir, probe, gres = sys.argv[1:5]
pr = pd.read_csv(probe, sep="\t").set_index("id")
g = pd.read_csv(gres, sep="\t").set_index("query_id")
rows = []
for i in pr.index:
    p = os.path.join(skdir, i + ".gsk")
    if not os.path.exists(p): continue
    c = pr.kmer_cov[i]
    a = read_abunds(p)
    sh, ex = features(a, c) if c >= MINC else (np.nan, np.nan)
    rows.append(dict(id=i, c=c, rust_cut=pr.rust_min_abund[i], est_mb=pr.est_mb[i], shoulder=sh, excess=ex))
d = pd.DataFrame(rows).set_index("id")
d["size_z"] = (d.est_mb - d.est_mb.median()) / (1.4826 * (d.est_mb - d.est_mb.median()).abs().median())
d["shoulder_z"] = robust_z(d.shoulder, d.c); d["excess_z"] = robust_z(d.excess, d.c)
d["second_frac"] = g.second_explained_fraction.reindex(d.index).fillna(0.0)
d["second_taxon"] = g.second_taxon.reindex(d.index)
d.to_csv(f"spectrum_{ds}.tsv", sep="\t")
ass = d[d.c >= MINC]
print(f"\n== {ds}: {len(d)} samples; assessable for the spectrum (c >= {MINC:g}): {len(ass)}; "
      f"index-confirmed mixtures: {int((d.second_frac > 0).sum())} (>=0.05: {int((d.second_frac >= 0.05).sum())}); "
      f"among assessable: {int((ass.second_frac > 0).sum())} (>=0.05: {int((ass.second_frac >= 0.05).sum())})")
for thr in (0.01, 0.05):
    y = (ass.second_frac >= thr).astype(int)
    if y.sum() < 2:
        print(f"  second_frac >= {thr}: only {int(y.sum())} assessable positive(s), AUC not computed"); continue
    out = []
    for col in ("shoulder_z", "excess_z", "size_z"):
        v = ass[col].fillna(ass[col].min() - 1) if col != "size_z" else ass[col]
        out.append(f"{col} AUC {roc_auc_score(y, v):.3f}")
    print(f"  positives = second_frac >= {thr} (n={int(y.sum())} of {len(y)}): " + "; ".join(out))
print("  assessable mixtures, with the spectrum statistics (z; cohort-relative):")
print(ass[ass.second_frac > 0].sort_values("second_frac", ascending=False)[["c", "second_taxon", "second_frac", "shoulder_z", "excess_z", "size_z"]].round(2).to_string())
for col in ("shoulder_z", "excess_z"):
    neg = ass[ass.second_frac == 0][col].dropna()
    print(f"  {col}: negatives above z=3: {int((neg > 3).sum())} of {len(neg)}; above z=5: {int((neg > 5).sum())}; "
          f"mixtures above z=3: {int((ass[ass.second_frac > 0][col] > 3).sum())} of {int((ass.second_frac > 0).sum())}")
nonass = d[(d.c < MINC) & (d.second_frac > 0)]
print(f"  mixtures NOT assessable by the spectrum (c < {MINC:g}): {len(nonass)}; their size_z:"); print(nonass[["c", "second_taxon", "second_frac", "size_z"]].round(2).to_string())
