# GmG against skani, Mash and sourmash (#46)

**Question.** Does a small FracMinHash index with per-rank voting beat the obvious alternative, "take the nearest reference genome and read its lineage", and what does it cost? Kraken2 and an ITS classifier were left out by decision (#46); this compares the four genome-level tools you would otherwise reach for.

**Setup.** Every tool searches the *same* reference: the ~7.4k genomes of the reps set of release 2026.10.2 (decoys removed), so the comparison is like-for-like and not full-BFD against reps. GmG runs `gmg run` with that reps index (scaled 20000). skani (`search`), Mash (`dist`, s=10000, k=21), Mash `screen` (`-w -i 0.75`) and sourmash (`gather`, k=21, scaled 10000) are each scored two ways with one code path (`scripts/compare/score.py`):

- *no cutoff*: the lineage of the best hit at every rank (what "nearest genome" gives);
- *GmG ANI cutoffs*: a rank is called only when the best hit's ANI reaches GmG's per-rank cutoff (0.7109 for phylum to order, 0.7715 family, 0.90 genus, 0.9938 species).

Leave-out runs drop the hits the query may not see *after* the search (itself, its species, its genus), which is valid for nearest-genome tools. Query sets: s1, 300 post-freeze NCBI assemblies (independent of the reference); s2, 300 RefSeq truth-set genomes scored leave-one/species/genus-out; s3, stress inputs (contamination 5/20/50%, 5-kb fragments, negative controls). Numbers are **precision of the calls made (share of inputs called)**. Each tool's ANI is its own estimate (skani ANI, 1 - Mash distance, sourmash containment ANI), comparable in size but not identical. Raw numbers: `tool_comparison.tsv`; reproduce with `scripts/compare/`.

## Tables


**Table 1. Post-freeze assemblies (s1, n = 300; independent of the reference). Precision of the calls made (share of inputs called).**

| tool | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|
| GmG (own calls) | 0.982 (0.95) | 0.985 (0.94) | 0.978 (0.93) | 0.967 (0.92) | 0.975 (0.83) | 1.000 (0.37) |
| skani nearest genome, no cutoff | 1.000 (0.83) | 1.000 (0.83) | 0.975 (0.83) | 0.980 (0.82) | 0.978 (0.80) | 0.654 (1.00) |
| skani nearest genome, GmG ANI cutoffs | 1.000 (0.83) | 1.000 (0.83) | 0.975 (0.83) | 0.980 (0.82) | 0.977 (0.78) | 0.778 (0.52) |
| mash dist nearest genome, no cutoff | 0.989 (0.91) | 0.989 (0.91) | 0.981 (0.91) | 0.967 (0.91) | 0.936 (0.89) | 0.500 (1.00) |
| mash dist nearest genome, GmG ANI cutoffs | 0.989 (0.91) | 0.989 (0.91) | 0.981 (0.91) | 0.970 (0.90) | 0.978 (0.81) | 0.741 (0.52) |
| mash screen nearest genome, no cutoff | 0.983 (0.95) | 0.967 (0.94) | 0.911 (0.94) | 0.897 (0.95) | 0.871 (0.93) | 0.500 (1.00) |
| mash screen nearest genome, GmG ANI cutoffs | 0.983 (0.95) | 0.967 (0.94) | 0.911 (0.94) | 0.901 (0.95) | 0.978 (0.80) | 0.741 (0.52) |
| sourmash gather nearest genome, no cutoff | 0.940 (0.78) | 0.969 (0.79) | 0.960 (0.78) | 0.904 (0.77) | 0.888 (0.79) | 0.765 (0.98) |
| sourmash gather nearest genome, GmG ANI cutoffs | 0.940 (0.78) | 0.969 (0.79) | 0.960 (0.78) | 0.920 (0.75) | 0.983 (0.40) | 0.944 (0.35) |


**Table 2a. RefSeq truth-set genomes, leave-one-out (s2). Precision (share called).**

| tool | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|
| GmG (own calls) | 0.960 (0.90) | 0.954 (0.92) | 0.949 (0.91) | 0.964 (0.86) | 0.988 (0.60) | 0.885 (0.19) |
| skani nearest genome, no cutoff | 1.000 (0.66) | 1.000 (0.67) | 1.000 (0.67) | 1.000 (0.67) | 0.972 (0.65) | 0.453 (0.66) |
| skani nearest genome, GmG ANI cutoffs | 1.000 (0.66) | 1.000 (0.67) | 1.000 (0.67) | 1.000 (0.67) | 0.981 (0.57) | 0.846 (0.19) |
| mash dist nearest genome, no cutoff | 0.987 (0.86) | 0.982 (0.87) | 0.978 (0.87) | 0.973 (0.87) | 0.875 (0.85) | 0.339 (0.86) |
| mash dist nearest genome, GmG ANI cutoffs | 0.987 (0.86) | 0.982 (0.87) | 0.978 (0.87) | 0.990 (0.80) | 0.982 (0.59) | 0.833 (0.22) |
| mash screen nearest genome, no cutoff | 0.886 (0.93) | 0.646 (0.95) | 0.543 (0.95) | 0.524 (0.96) | 0.502 (0.92) | 0.302 (0.93) |
| mash screen nearest genome, GmG ANI cutoffs | 0.886 (0.93) | 0.646 (0.95) | 0.543 (0.95) | 0.525 (0.94) | 0.989 (0.34) | 0.942 (0.19) |
| sourmash gather nearest genome, no cutoff | 0.936 (0.80) | 0.925 (0.83) | 0.925 (0.83) | 0.886 (0.82) | 0.788 (0.79) | 0.173 (0.80) |
| sourmash gather nearest genome, GmG ANI cutoffs | 0.936 (0.80) | 0.925 (0.83) | 0.925 (0.83) | 0.920 (0.78) | 0.975 (0.44) | 0.714 (0.10) |


**Table 2b. RefSeq truth-set genomes, leave-species-out (s2). Precision (share called).**

| tool | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|
| GmG (own calls) | 0.953 (0.86) | 0.948 (0.90) | 0.943 (0.88) | 0.957 (0.81) | 0.975 (0.44) | 0.000 (0.03) |
| skani nearest genome, no cutoff | 1.000 (0.54) | 1.000 (0.56) | 1.000 (0.56) | 1.000 (0.56) | 0.952 (0.53) | 0.000 (0.55) |
| skani nearest genome, GmG ANI cutoffs | 1.000 (0.54) | 1.000 (0.56) | 1.000 (0.56) | 1.000 (0.56) | 0.964 (0.41) | 0.000 (0.04) |
| mash dist nearest genome, no cutoff | 0.987 (0.84) | 0.982 (0.86) | 0.977 (0.86) | 0.973 (0.86) | 0.846 (0.83) | 0.000 (0.84) |
| mash dist nearest genome, GmG ANI cutoffs | 0.987 (0.84) | 0.982 (0.86) | 0.977 (0.86) | 0.990 (0.75) | 0.966 (0.43) | 0.000 (0.05) |
| mash screen nearest genome, no cutoff | 0.832 (0.91) | 0.512 (0.95) | 0.391 (0.95) | 0.365 (0.94) | 0.320 (0.89) | 0.000 (0.92) |
| mash screen nearest genome, GmG ANI cutoffs | 0.832 (0.91) | 0.512 (0.95) | 0.391 (0.95) | 0.360 (0.92) | 0.947 (0.07) | 0.000 (0.01) |
| sourmash gather nearest genome, no cutoff | 0.923 (0.71) | 0.911 (0.74) | 0.911 (0.75) | 0.868 (0.74) | 0.751 (0.70) | 0.000 (0.71) |
| sourmash gather nearest genome, GmG ANI cutoffs | 0.923 (0.71) | 0.911 (0.74) | 0.911 (0.75) | 0.904 (0.69) | 0.968 (0.34) | 0.000 (0.03) |


**Table 2c. RefSeq truth-set genomes, leave-genus-out (s2). Precision (share called).**

| tool | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|
| GmG (own calls) | 0.820 (0.71) | 0.776 (0.74) | 0.681 (0.65) | 0.615 (0.51) | 0.000 (0.07) | - (0.00) |
| skani nearest genome, no cutoff | 1.000 (0.15) | 0.976 (0.16) | 0.714 (0.16) | 0.700 (0.16) | 0.000 (0.12) | 0.000 (0.16) |
| skani nearest genome, GmG ANI cutoffs | 1.000 (0.15) | 0.976 (0.16) | 0.714 (0.16) | 0.700 (0.16) | 0.000 (0.05) | 0.000 (0.01) |
| mash dist nearest genome, no cutoff | 0.933 (0.49) | 0.895 (0.52) | 0.771 (0.51) | 0.713 (0.47) | 0.000 (0.42) | 0.000 (0.49) |
| mash dist nearest genome, GmG ANI cutoffs | 0.933 (0.49) | 0.895 (0.52) | 0.771 (0.51) | 0.739 (0.36) | 0.000 (0.06) | 0.000 (0.01) |
| mash screen nearest genome, no cutoff | 0.759 (0.88) | 0.318 (0.90) | 0.112 (0.90) | 0.056 (0.90) | 0.000 (0.86) | 0.000 (0.88) |
| mash screen nearest genome, GmG ANI cutoffs | 0.759 (0.88) | 0.318 (0.90) | 0.112 (0.90) | 0.054 (0.86) | 0.000 (0.00) | 0.000 (0.01) |
| sourmash gather nearest genome, no cutoff | 0.667 (0.38) | 0.606 (0.40) | 0.558 (0.40) | 0.410 (0.39) | 0.000 (0.37) | 0.000 (0.39) |
| sourmash gather nearest genome, GmG ANI cutoffs | 0.667 (0.38) | 0.606 (0.40) | 0.558 (0.40) | 0.448 (0.34) | 0.000 (0.04) | - (0.00) |


**Table 3. Stress inputs (s3): genus precision (share called) per degradation; the last column is the share of negative controls that got a phylum call (lower is better).**

| tool | complete | contam:20 | contam:5 | contam:50 | frag5k | negatives called |
|---|---|---|---|---|---|---|
| GmG (own calls) | 1.000 (0.79) | 0.792 (1.00) | 1.000 (0.79) | 0.773 (0.92) | 1.000 (0.79) | 0.00 |
| skani nearest genome, no cutoff | 1.000 (0.88) | 0.708 (1.00) | 1.000 (0.88) | 0.667 (1.00) | 1.000 (0.88) | 0.00 |
| skani nearest genome, GmG ANI cutoffs | 1.000 (0.79) | 0.708 (1.00) | 1.000 (0.79) | 0.667 (1.00) | 1.000 (0.83) | 0.00 |
| mash dist nearest genome, no cutoff | 1.000 (0.95) | 0.792 (1.00) | 0.875 (1.00) | 0.750 (1.00) | 1.000 (0.96) | 0.00 |
| mash dist nearest genome, GmG ANI cutoffs | 1.000 (0.79) | 0.792 (1.00) | 1.000 (0.79) | 0.750 (1.00) | 1.000 (0.79) | 0.00 |
| mash screen nearest genome, no cutoff | 0.947 (0.99) | 0.792 (1.00) | 0.875 (1.00) | 0.667 (1.00) | 0.875 (1.00) | 0.00 |
| mash screen nearest genome, GmG ANI cutoffs | 1.000 (0.61) | 0.792 (1.00) | 1.000 (0.79) | 0.667 (1.00) | 1.000 (0.79) | 0.00 |
| sourmash gather nearest genome, no cutoff | 0.904 (0.86) | 0.739 (0.96) | 0.783 (0.96) | 0.609 (0.96) | 0.905 (0.88) | 0.00 |
| sourmash gather nearest genome, GmG ANI cutoffs | 0.941 (0.71) | 0.789 (0.79) | 0.941 (0.71) | 0.684 (0.79) | 0.941 (0.71) | 0.00 |


**Table 4. Cost (8 CPUs per run; wall seconds per input over the whole run, including loading the reference once; peak memory).**

| tool | s1 s/input | s2 s/input | s3 s/input | peak RAM GB | reference build (wall) | reference size |
|---|---|---|---|---|---|---|
| GmG | 0.13 | 0.13 | 0.18 | 0.4 | (index build ~10 min on the full BFD; reps index from the release) | 67 MB |
| skani | 1.34 | 0.80 | 1.75 | 23.9 | 17 min | 30740 MB |
| mash dist | 0.71 | 0.70 | 0.78 | 1.2 | 25 min | 594 MB |
| mash screen | 33.74 | 30.93 | 30.38 | 13.3 | 25 min | 594 MB |
| sourmash gather | 5.23 | 4.85 | 4.88 | 0.3 | 248 min | 89 MB |

## What the tables say

Read precision together with the share called: a tool that abstains on hard inputs looks precise. "Correct per input" is precision x called share.

1. **Known genomes (s1, s2 leave-one-out): the tools are close, nothing clearly dominates.** At genus GmG is 0.975 (0.83 called) on post-freeze assemblies and 0.988 (0.60) on leave-one-out truth genomes; skani with GmG cutoffs is 0.977 (0.78) and 0.981 (0.57), Mash `dist` 0.978 (0.81) and 0.982 (0.59). At phylum to family Mash `dist` and skani are as precise as GmG or a little better (leave-one-out phylum 0.987 for Mash dist and 1.000 for skani against GmG's 0.960), but skani only calls 0.66 of those inputs (it reports no hit when nothing is near) where GmG calls 0.90. So on genomes that have a relative in the reference, a 67 MB sketch index loses nothing material at genus and costs a few points of precision at the high ranks in exchange for calling more inputs.
2. **Novel genera (s2 leave-genus-out): this is where GmG's voting pays.** Correct class calls per input: GmG 0.58 (0.776 x 0.74), Mash `dist` 0.47 (0.895 x 0.52), Mash `screen` 0.29, sourmash 0.24, skani 0.16 (precision 0.976 but only 0.16 called). Phylum: GmG 0.58, Mash `dist` 0.46, skani 0.15. GmG's precision is lower than Mash dist's (0.776 against 0.895 at class), so the gain is coverage, not accuracy: it keeps calling the higher ranks where the nearest-genome tools find no hit above their own threshold. At order and family all tools degrade (GmG 0.68 and 0.62); no tool should be trusted for a new genus below class.
3. **Species: GmG abstains where the others guess.** Without a cutoff, every nearest-genome tool names a species for every input and is right only part of the time (post-freeze 0.50-0.77; leave-one-out 0.17-0.45). With GmG's 0.9938 cutoff Mash and skani reach 0.74-0.78 at 0.52 called on post-freeze assemblies and sourmash 0.94 at 0.35; GmG gets 1.000 at 0.37 called, with no wrong species among the inputs it named. The cost is a lower called share (0.37), by design.
4. **Contamination (s3).** Genus precision on 20% and 50% contaminated inputs is 0.79 and 0.77 for GmG (0.92 called at 50%), 0.79 and 0.75 for Mash `dist`, 0.71 and 0.67 for skani, 0.74-0.79 and 0.61-0.68 for sourmash: GmG is no worse than any of them and ahead of skani and sourmash, with a lower share called than the no-cutoff tools. No tool called a phylum for any negative control, but none had decoys in its reference, so that column only shows that the negatives are far from every fungal genome; it does not test non-fungal rejection (see `decoys.md`).
5. **Cost: this is the clear difference.**


| | GmG | Mash `dist` | skani | sourmash `gather` | Mash `screen` |
|---|---|---|---|---|---|
| seconds per input (8 CPUs, incl. loading) | 0.13-0.18 | 0.7 | 0.8-1.8 | 4.9-5.2 | 30-34 |
| peak RAM | 0.4 GB | 1.2 GB | 23.9 GB | 0.3 GB | 13.3 GB |
| reference on disk | 67 MB | 594 MB | 30.7 GB | 89 MB | 594 MB |
| reference build | prebuilt (release asset) | 25 min | 17 min | 248 min | 25 min |

GmG is 5x faster than Mash `dist`, 6-13x faster than skani, ~35x faster than sourmash and ~200x faster than Mash `screen`, with a reference 9-450x smaller than Mash's or skani's. skani's 24 GB RAM and 31 GB database are what rule it out for a laptop or the browser.

## Caveats

- Same ~7.4k-genome reference for everyone. A larger reference (the full 23.7k BFD) would help every nearest-genome tool, GmG's full index included, and was not run here.
- "Nearest genome" ranking uses each tool's own ANI; the GmG cutoffs were calibrated on GmG's containment ANI, so the cutoff variants are generous to the others at some ranks and strict at others (see the species rows).
- s1 and s2 are 300 genomes each; one input is 0.3 points of share, and the stress groups have ~24 inputs (one input is ~4 points). Differences under about 3 points are noise.
- Leave-genus-out is applied after the search, which is exact for the nearest-genome tools. It is also what the GmG `--exclude-level genus` run does (`reproducibility.md`).
- Mash `screen` is shown for completeness; its containment identity is not designed to rank genomes, which is why its no-cutoff precision falls fast below class.

## Conclusion

GmG is not more accurate than a good nearest-genome search on genomes that have a close relative in the reference: Mash `dist` and skani are within a few points at genus and sometimes ahead at phylum. It earns its place in two ways. It calls more inputs correctly when the genus is new (class 0.58 correct per input against 0.47 for the best alternative, 0.16 for skani), and it is far cheaper (0.4 GB, 67 MB, 0.13 s per input against 24 GB, 31 GB and 1.3 s for skani; Mash `dist` is the nearest on cost). It also abstains at species where uncut nearest-genome tools are right only 50-77% of the time. Nothing here changes the defaults. Kraken2 and an ITS classifier remain out of scope (#46); the full-BFD comparison is the obvious next run if the question matters for the paper.
