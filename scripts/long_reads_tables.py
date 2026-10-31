#!/usr/bin/env python3
"""Markdown pivot tables (error rate x depth) from the accuracy table of long_reads_score.py summarize (#44).

usage: long_reads_tables.py ACCURACY.tsv > tables.md
"""
import csv, sys

A = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
ERRS = [0.5, 1, 2, 5, 8, 10, 12, 15]
DEPS = [0.25, 0.5, 1, 2, 5, 10, 30]


def table(title, tag, rank, col, fmt="{:.2f}"):
    print(f"**{title}**\n")
    print("| error \\ depth | " + " | ".join(f"{d:g}x" for d in DEPS) + " |")
    print("|---|" + "---|" * len(DEPS))
    for e in ERRS:
        cells = []
        for d in DEPS:
            r = [x for x in A if x["tag"] == tag and x["rank"] == rank and float(x["nominal_error_pct"]) == e and float(x["depth"]) == d]
            try:
                cells.append(fmt.format(float(r[0][col])))
            except (IndexError, ValueError):
                cells.append("-")
        print(f"| {e:g}% | " + " | ".join(cells) + " |")
    print()


T = [
    ("Genus called (share of the 24 genomes), own genome hidden, min-abund 2", "self_a2", "genus", "called"),
    ("Genus precision of the calls made, own genome hidden, min-abund 2", "self_a2", "genus", "precision"),
    ("Phylum called, own genome hidden, min-abund 2", "self_a2", "phylum", "called"),
    ("Species recall over the 17 genomes with a conspecific in the reference, min-abund 2 (default)", "self_a2", "species", "recall"),
    ("Species recall, min-abund 1", "self_a1", "species", "recall"),
    ("Species recall, min-abund auto", "self_aauto", "species", "recall"),
    ("Species precision of the calls made, min-abund 2", "self_a2", "species", "precision"),
    ("Genus called with the species hidden (genus precision follows), min-abund 2", "species_a2", "genus", "called"),
    ("Genus precision with the species hidden", "species_a2", "genus", "precision"),
    ("Genus called with the genus hidden (every call is wrong), min-abund 2", "genus_a2", "genus", "called"),
    ("Median ANI of the top hit, own genome hidden, min-abund 2", "self_a2", "species", "median_ani"),
    ("Median ANI, min-abund 1", "self_a1", "species", "median_ani"),
    ("Median `coverage_est` over the true depth, min-abund 2", "self_a2", "species", "median_cov_est_over_depth"),
]
for t in T:
    table(*t, fmt="{:.3f}" if "ani" in t[3] else "{:.2f}")
