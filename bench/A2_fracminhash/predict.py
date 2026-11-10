#!/usr/bin/env python3
"""A2, prediction: sampling noise and saturation of the FracMinHash containment ANI. Arithmetic only (no data).

Model. A genome has about G distinct k-mers; a FracMinHash sketch at `scaled` s keeps those with hash <= 2^64 / s, so M = G / s
hashes. Two genomes share a fraction C = ANI^k of their k-mers. The estimate is C_hat = X / M with X ~ Binomial(M, C), and
ANI_hat = C_hat^(1/k), so  SD(ANI_hat) = (ANI / k) * sqrt((1 - C) / (C * M))  (delta method).
Saturation: the estimate needs shared hashes; with X = M * C below about MIN_SHARED the estimate is unusable.
usage: predict.py   (prints tables as TSV)"""
import math

K = 21
MIN_SHARED = 10  # shared sketch hashes below which the estimate is not trusted


def sd_ani(ani, g, s, k=K):
    c = ani**k
    m = g / s
    return (ani / k) * math.sqrt((1 - c) / (c * m))


def lowest_ani(g, s, k=K, min_shared=MIN_SHARED):
    """ANI at which the expected number of shared hashes drops to min_shared"""
    m = g / s
    return (min_shared / m) ** (1 / k)


def main():
    print("table\tgenome_mb\tscaled\tani\tsd_1_minus_ani\thashes")
    for g in (20e6, 40e6, 140e6):
        for s in (200, 1000):
            for ani in (0.9999, 0.999, 0.99, 0.95, 0.9, 0.8, 0.7):
                print(f"noise\t{g/1e6:.0f}\t{s}\t{ani}\t{sd_ani(ani, g, s):.2e}\t{g/s:.0f}")
    print("table\tgenome_mb\tscaled\tlowest_ani_with_10_shared\tlowest_ani_with_30_shared")
    for g in (20e6, 40e6, 140e6):
        for s in (200, 1000):
            print(f"saturation\t{g/1e6:.0f}\t{s}\t{lowest_ani(g, s):.3f}\t{lowest_ani(g, s, min_shared=30):.3f}")


if __name__ == "__main__":
    main()
