#!/usr/bin/env python3
"""26_gap_null.py : distribution of the largest-gap ratio (gap1 / gap2) for a population with no structure.

For a Kingman coalescent the UPGMA merge heights of the true tree are cumulative sums of independent waiting times
T_k ~ Exp(k(k-1)/2) (k lineages), so the gap at the cut that leaves k clusters is T_k. The ratio is scale-free. This gives the
chance of calling a 'clear cut' (ratio above a threshold) in a structureless population, for a given number of samples."""
import numpy as np
rng = np.random.default_rng(1)
def null_ratios(n, kmax=30, reps=20000):
    ks = np.arange(2, min(kmax, n) + 1)
    rate = ks * (ks - 1) / 2.0
    T = rng.exponential(1.0 / rate, size=(reps, len(ks)))
    top2 = np.sort(T, axis=1)[:, -2:]
    return top2[:, 1] / top2[:, 0]
print("null distribution of gap1/gap2 for a structureless (Kingman coalescent) population, kmax=30")
print("  n   median   90th    95th    99th   P(ratio>=2)  P(>=3)  P(>=5)  P(>=10)")
for n in (10, 30, 100, 400):
    r = null_ratios(n)
    print(f"{n:4d}  {np.median(r):6.2f}  {np.quantile(r,.9):6.2f}  {np.quantile(r,.95):6.2f}  {np.quantile(r,.99):6.2f}   "
          f"{(r>=2).mean():6.3f}  {(r>=3).mean():6.3f}  {(r>=5).mean():6.3f}  {(r>=10).mean():6.3f}")
print("\nobserved (Rust cutoff reruns): Cocci both species 8.9, immitis 69.2 (one broken sample), posadasii 6.0, Bd 1.96")
