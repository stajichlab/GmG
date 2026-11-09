#!/usr/bin/env python3
"""Figures for reports/fusarium_markers/README.md, from the result tables committed in the repo.

Reads bench/F1_fusariumid/results/{markers_multilocus,markers_nearest_species,markers_kmer,genomes_bench_summary}.tsv and
bench/F2_marker_mode/results/calibration.tsv; writes PNGs into reports/fusarium_markers/figures/.
Palette: the dataviz reference categorical slots 1-4 (validated light mode: CVD and normal-vision separation pass;
slots 3-4 are below 3:1 against the surface, so every series is also direct-labelled and has its own marker/line style).
usage: python3 reports/fusarium_markers/make_figures.py   (run from the repository root)
"""
import collections
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
F1 = os.path.join(ROOT, "bench", "F1_fusariumid", "results")
F2 = os.path.join(ROOT, "bench", "F2_marker_mode", "results")
OUT = os.path.join(ROOT, "reports", "fusarium_markers", "figures")
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # categorical slots 1-4, fixed order
MARK = ["o", "s", "^", "D"]
DASH = ["-", "--", "-.", ":"]


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def style(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURFACE)
    ax.figure.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_xlabel(xlabel, color=INK2, fontsize=10)
    ax.set_ylabel(ylabel, color=INK2, fontsize=10)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=150, facecolor=SURFACE)
    plt.close(fig)


def loci_accuracy():
    r = rows(os.path.join(F1, "markers_multilocus.tsv"))
    r = sorted(r, key=lambda x: float(x["accuracy"]))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    y = range(len(r))
    acc = [float(x["accuracy"]) for x in r]
    ax.barh(list(y), acc, color=S[0], height=0.6, edgecolor=SURFACE, linewidth=2)
    for i, x in zip(y, r):
        ax.text(float(x["accuracy"]) + 0.01, i, f"{float(x['accuracy']):.0%}  (n={int(x['answerable']):,})",
                va="center", color=INK2, fontsize=8.5)
    ax.set_yticks(list(y), [x["set"].replace("+", " + ") for x in r], color=INK, fontsize=9)
    ax.set_xlim(0, 1.18)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0], ["0%", "25%", "50%", "75%", "100%"])
    style(ax, "Leave-one-out species accuracy by marker set",
          "nearest other strain is the same species (answerable queries)", "")
    ax.grid(axis="y", visible=False)
    save(fig, "loci_accuracy.png")


def nearest_diffs():
    r = rows(os.path.join(F1, "markers_nearest_species.tsv"))
    sets = ["tef1", "tef1+rpb2", "tef1+rpb1+rpb2", "tef1+rpb1+rpb2+CaM+tub2"]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for i, s in enumerate(sets):
        d = sorted(int(x["min_diffs_to_other_species"]) for x in r if x["set"] == s)
        xs = list(range(0, 61))
        ys = [sum(v <= t for v in d) / len(d) for t in xs]
        ax.step(xs, ys, where="post", color=S[i], linewidth=2, linestyle=DASH[i], label=f"{s} ({len(d)} species)")
        ax.annotate(s.replace("+", " + "), (60, ys[-1]), xytext=(4, 0), textcoords="offset points", va="center",
                    color=INK2, fontsize=8.5)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0], ["0%", "25%", "50%", "75%", "100%"])
    style(ax, "How little sequence separates a species from its closest relative",
          "differences (mismatches + gap openings) to the nearest strain of another species",
          "species within that many differences")
    ax.legend(frameon=False, fontsize=8.5, loc="lower right", labelcolor=INK2)
    save(fig, "nearest_species_diffs.png")


def kmer_scaled():
    r = rows(os.path.join(F1, "markers_kmer.tsv"))
    loci = ["tef1", "rpb2", "rpb1", "ITS"]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for i, l in enumerate(loci):
        pts = sorted((float(x["mean_hashes_per_seq"]), float(x["accuracy"])) for x in r if x["locus"] == l)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=S[i], linewidth=2, linestyle=DASH[i],
                marker=MARK[i], markersize=8, markeredgecolor=SURFACE, markeredgewidth=2, label=l)
        ax.annotate(l, pts[-1], xytext=(6, 0), textcoords="offset points", va="center", color=INK2, fontsize=9)
    ax.set_xscale("log")
    ax.set_ylim(0.2, 1)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0], ["20%", "40%", "60%", "80%", "100%"])
    ax.axvline(0.64, color=INK2, linewidth=1, linestyle=":")
    ax.annotate("tef1 at genome-index\nscaled 1000: ~0.6 hashes", (0.64, 0.3), xytext=(6, 0), textcoords="offset points",
                color=INK2, fontsize=8)
    ax.set_xlim(0.4, 2500)
    style(ax, "k-mer nearest neighbour on one marker: how many hashes are needed",
          "k-mer hashes kept per sequence (FracMinHash scaled 1, 3, 10, 30, 100; log scale)", "species accuracy (LOO)")
    ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK2)
    save(fig, "kmer_hashes.png")


def marker_mode_cutoff():
    r = [x for x in rows(os.path.join(F2, "calibration.tsv")) if x["locus"] == "all"]
    c = [float(x["cutoff"]) for x in r]
    series = [("species precision", "species_precision"), ("species called", "species_called"),
              ("absent species named (false call)", "absent_false_call")]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for i, (lab, k) in enumerate(series):
        ys = [float(x[k]) for x in r]
        ax.plot(c, ys, color=S[i], linewidth=2, linestyle=DASH[i], marker=MARK[i], markersize=8,
                markeredgecolor=SURFACE, markeredgewidth=2, label=lab)
        ax.annotate(lab, (c[-1], ys[-1]), xytext=(6, 0), textcoords="offset points", va="center", color=INK2,
                    fontsize=8.5)
    ax.axvline(0.995, color=INK2, linewidth=1, linestyle=":")
    ax.annotate("shipped: 0.995", (0.995, 0.9), xytext=(-70, 0), textcoords="offset points", color=INK2, fontsize=8.5)
    ax.set_ylim(0, 1.02)
    ax.set_xlim(0.979, 1.0085)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0], ["0%", "25%", "50%", "75%", "100%"])
    style(ax, "Marker mode: choosing the species ANI cutoff (3,408 held-out marker queries)",
          "species ANI cutoff", "fraction of queries")
    ax.legend(frameon=False, fontsize=8.5, loc="center left", labelcolor=INK2)
    save(fig, "marker_mode_cutoff.png")


def genome_vs_markers():
    g = [x for x in rows(os.path.join(F1, "genomes_bench_summary.tsv")) if x["refset"] == "fus" and x["mode"] == "loo"]
    g = sorted(g, key=lambda x: int(x["median_query_hashes"]))
    k = sorted((float(x["mean_hashes_per_seq"]), float(x["accuracy"])) for x in rows(os.path.join(F1, "markers_kmer.tsv"))
               if x["locus"] == "tef1")
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    gx = [int(x["median_query_hashes"]) for x in g]
    gy = [float(x["fid_top_correct"]) for x in g]
    ax.plot(gx, gy, color=S[0], linewidth=2, marker=MARK[0], markersize=8, markeredgecolor=SURFACE, markeredgewidth=2,
            label="whole genome (2,110 BFD genomes)")
    ax.annotate("whole genome", (gx[-1], gy[-1]), xytext=(6, 0), textcoords="offset points", va="center", color=INK2,
                fontsize=9)
    ax.plot([p[0] for p in k], [p[1] for p in k], color=S[1], linewidth=2, linestyle="--", marker=MARK[1], markersize=8,
            markeredgecolor=SURFACE, markeredgewidth=2, label="tef1 marker alone (FusariumID)")
    ax.annotate("tef1 alone", k[-1], xytext=(6, 0), textcoords="offset points", va="center", color=INK2, fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(4, 2e5)
    ax.set_ylim(0.6, 1)
    ax.set_yticks([0.6, 0.7, 0.8, 0.9, 1.0], ["60%", "70%", "80%", "90%", "100%"])
    style(ax, "Species accuracy vs hashes per query: genomes have information to spare",
          "k-mer hashes per query (log scale)", "nearest reference is the same FusariumID species")
    ax.legend(frameon=False, fontsize=9, loc="lower right", labelcolor=INK2)
    save(fig, "genome_vs_marker_hashes.png")


def main():
    os.makedirs(OUT, exist_ok=True)
    loci_accuracy()
    nearest_diffs()
    kmer_scaled()
    marker_mode_cutoff()
    genome_vs_markers()
    print(f"figures -> {OUT}: {sorted(os.listdir(OUT))}")


if __name__ == "__main__":
    main()
