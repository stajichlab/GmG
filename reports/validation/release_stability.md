# Stability under reference and taxonomy updates (#48)

Four releases, classified by the **same binary** (`gmg-a9acf5c`, `known` thresholds) so that only the index differs:

| release | what changed against the one before |
|---|---|
| 2026.10 | BFD, 23,669 genomes, BFD species strings as labels |
| 2026.10.1 | NCBI taxids and a normalised NCBI taxonomy (taxdump 2026-10-29) replace the BFD strings; the same genomes |
| 2026.10.2 | + 81 non-fungal decoy genomes (#79) |
| 2026.10.3 | - 57 genomes with ambiguous labels (`reports/label_audit/reference_exclusions.txt`); decoys kept |

Test sets: the 1,257 post-freeze NCBI assemblies (truth = NCBI lineage; `data/validation/post_freeze.meta.tsv`), the 199 IndigoAg assemblies (no lineage truth: changes only),
the 15 held-out non-fungal genomes. `scripts/hpcc/25_release_stability.sbatch` and `26_release_stability_small.sbatch` run them, `scripts/release_compare.py` compares them and applies the
re-calibration rules of `docs/recalibration.md`. Per-genome results stay on the HPCC (`GmG_runs/external_big/rel_*`, `GmG_runs/relstab`).

## Findings

* **Updates move few calls, and the cause of each move is known.** Of the 1,257 post-freeze assemblies: 2026.10 to .1 changed **8 species labels (0.6%)** and nothing else (taxonomy
  normalisation: *Zygosaccharomyces sapae* to *Z. rouxii*, *Emmonsia* bracket labels; the same top hit, a new label); .1 to .2 changed 1 call (0.1%, a label); .2 to .3 changed 3 (0.2%: one lost its excluded
  top hit, two label updates; family and genus calls). **No genus call among the 1,045 confident (`ok`) inputs changed in any step.**
* **The calibration holds.** Precision of the calls made on the post-freeze set is the same in all four releases to three decimals: phylum 0.990, class 0.991, order 0.986, family 0.983-0.984, genus 0.988,
  species 1.000; the called share is the same (phylum 0.95, genus 0.84, species 0.75; family 0.90 to 0.91). The decoys (2026.10.2) cost nothing on fungal calls *because of the distant-decoy rule* (#89): the first decoy build,
  before it, put a Chordata phylum on 25 post-freeze fungi and lowered phylum precision from 0.990 to 0.977 (`reports/validation/decoys.md`, `truth_set.md`).
* **`weak_phylum` (#99) is stable across releases.** 70 of the 1,189 phylum calls on the post-freeze set (5.9%) come from a best-hit ANI below 0.85; phylum precision is **0.996 without them and 0.886 with them**, and the
  numbers are identical in all four releases (derived from `ani`, so the older runs are comparable). The flag keeps the call; `release_compare.py` now reports the share and the split, and a move of the share by more than 5 points is a
  re-calibration rule (rule 2 of `docs/recalibration.md`).
* **Decoys change what they should.** On the 15 held-out non-fungal genomes, 4 changed (the three that are themselves decoys and *Entamoeba*, which has a decoy relative: all now `not_fungal`) and 11 did not (bacteria, viruses and the non-fungi without a close decoy keep their calls); nothing else in the other sets moved with the decoys
  (1 IndigoAg call, a label).
* **Removing 57 genomes changes the inputs whose nearest genome they were, by design.** IndigoAg 2026.10.2 to .3: 5 inputs changed (2.5%): three assemblies called *Pseudopestalotiopsis* from an excluded genome (*P. chinensis*, one of the
  14 suspects) are now *Pestalotiopsis* (the provider's name; ANI 0.982 to *P. fici*); one assembly with a 0.99-ANI hit to an excluded *Tetracladium* genome loses its confident call and gets a low-confidence call
  from a distant, chance hit (ANI 0.81); one low-confidence phylum flips. The loss of the *Tetracladium* call is the cost of removing a genome that was probably mislabelled but was also the only close relative; this is the
  "distant relatives give wrong low-confidence phylum calls" weakness already seen in IndigoAg (`external_indigo.md`).
* **Cluster ids are stable until genomes are removed.** Clusters are computed once and carried: ids and partition identical across 2026.10, .1, .2 (adjusted Rand index 1.0000, 23,669 of 23,669 genomes). 2026.10.3 keeps the ids of
  its 23,612 genomes, but 15 removed genomes were centroids: **3 clusters (23 genomes) now carry the id of a genome that is not in the index** (`scripts/reseed_clusters.py` fixes this and the release script now calls it; the published `index-2026.10.3` predates it and still has them; see `clusters_vs_skani.md`).

## Re-calibration trigger

`docs/recalibration.md` defines five rules; `release_compare.py` checks the first three automatically (exit code 3 when one fires). **No rule fired for any release step** (the table below). The one confident-call change
that crossed 2% of a set in a first version of the rule (IndigoAg .2 to .3: 3 of 137, 2.2%) was caused by the deliberate exclusion and is therefore not counted: rule 3 counts unexplained changes only.
Rule 4 (golden panel) is `reports/validation/reproducibility.md`; rule 5 (size and taxonomy renames) did not apply (the largest step, removing 57 of 23,750 genomes, is 0.2%).

## Tables (generated by `scripts/release_compare.py`)


### held

| from | to | inputs | status changed | any call changed | phylum | class | order | family | genus | species | `ok` in both, genus or species call changed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2026.10 | 2026.10.1 | 15 | 0 | 0 (0.0%) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2026.10.1 | 2026.10.2 | 15 | 4 | 4 (26.7%) | 4 | 4 | 4 | 4 | 4 | 3 | 0 |
| 2026.10.2 | 2026.10.3 | 15 | 0 | 0 (0.0%) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Causes of the changed inputs (first matching rule):

| cause | 2026.10→2026.10.1 | 2026.10.1→2026.10.2 | 2026.10.2→2026.10.3 |
|---|---|---|---|
| a decoy genome is (no longer) the top hit | 0 | 4 | 0 |

### indigo

| from | to | inputs | status changed | any call changed | phylum | class | order | family | genus | species | `ok` in both, genus or species call changed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2026.10 | 2026.10.1 | 199 | 0 | 0 (0.0%) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 2026.10.1 | 2026.10.2 | 199 | 0 | 1 (0.5%) | 0 | 0 | 1 | 1 | 0 | 0 | 0 |
| 2026.10.2 | 2026.10.3 | 199 | 1 | 5 (2.5%) | 2 | 2 | 2 | 2 | 4 | 0 | 3 |

Causes of the changed inputs (first matching rule):

| cause | 2026.10→2026.10.1 | 2026.10.1→2026.10.2 | 2026.10.2→2026.10.3 |
|---|---|---|---|
| old top hit excluded from the reference | 0 | 0 | 4 |
| same top hit, new label (taxonomy update) | 0 | 1 | 1 |

### post

| from | to | inputs | status changed | any call changed | phylum | class | order | family | genus | species | `ok` in both, genus or species call changed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2026.10 | 2026.10.1 | 1257 | 0 | 8 (0.6%) | 0 | 0 | 0 | 0 | 0 | 8 | 8 |
| 2026.10.1 | 2026.10.2 | 1257 | 0 | 1 (0.1%) | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| 2026.10.2 | 2026.10.3 | 1257 | 3 | 3 (0.2%) | 0 | 0 | 0 | 3 | 3 | 0 | 0 |

Causes of the changed inputs (first matching rule):

| cause | 2026.10→2026.10.1 | 2026.10.1→2026.10.2 | 2026.10.2→2026.10.3 |
|---|---|---|---|
| old top hit excluded from the reference | 0 | 0 | 1 |
| same top hit, new label (taxonomy update) | 8 | 1 | 2 |

Precision of the calls made (called share), against the truth table:

| release | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|
| 2026.10 | 0.990 (0.95) | 0.991 (0.94) | 0.986 (0.94) | 0.983 (0.90) | 0.988 (0.84) | 1.000 (0.75) |
| 2026.10.1 | 0.990 (0.95) | 0.991 (0.94) | 0.986 (0.94) | 0.983 (0.90) | 0.988 (0.84) | 1.000 (0.75) |
| 2026.10.2 | 0.990 (0.95) | 0.991 (0.94) | 0.986 (0.94) | 0.984 (0.90) | 0.988 (0.84) | 1.000 (0.75) |
| 2026.10.3 | 0.990 (0.95) | 0.991 (0.94) | 0.986 (0.94) | 0.984 (0.91) | 0.988 (0.84) | 1.000 (0.75) |

Phylum calls with a best-hit ANI below 0.85 (`weak_phylum`, kept as calls and flagged): share of the phylum calls, and phylum precision of the calls without and with the flag:

| release | weak share | precision without | precision with (weak only) |
|---|---|---|---|
| 2026.10 | 70 of 1189 (5.9%) | 0.996 (1119) | 0.886 (70) |
| 2026.10.1 | 70 of 1189 (5.9%) | 0.996 (1119) | 0.886 (70) |
| 2026.10.2 | 70 of 1189 (5.9%) | 0.996 (1119) | 0.886 (70) |
| 2026.10.3 | 70 of 1189 (5.9%) | 0.996 (1119) | 0.886 (70) |

### Re-calibration trigger

Churn of confident genus calls (rule 3 counts only the unexplained ones):

* held 2026.10->2026.10.1: 0 of 0 confident genus calls changed, 0 of them explained (a deliberately excluded top hit or a label update), 0 unexplained
* held 2026.10.1->2026.10.2: 0 of 0 confident genus calls changed, 0 of them explained (a deliberately excluded top hit or a label update), 0 unexplained
* held 2026.10.2->2026.10.3: 0 of 4 confident genus calls changed, 0 of them explained (a deliberately excluded top hit or a label update), 0 unexplained
* indigo 2026.10->2026.10.1: 0 of 138 confident genus calls changed, 0 of them explained (a deliberately excluded top hit or a label update), 0 unexplained
* indigo 2026.10.1->2026.10.2: 0 of 138 confident genus calls changed, 0 of them explained (a deliberately excluded top hit or a label update), 0 unexplained
* indigo 2026.10.2->2026.10.3: 3 of 137 confident genus calls changed, 3 of them explained (a deliberately excluded top hit or a label update), 0 unexplained
* post 2026.10->2026.10.1: 0 of 1045 confident genus calls changed, 0 of them explained (a deliberately excluded top hit or a label update), 0 unexplained
* post 2026.10.1->2026.10.2: 0 of 1045 confident genus calls changed, 0 of them explained (a deliberately excluded top hit or a label update), 0 unexplained
* post 2026.10.2->2026.10.3: 0 of 1045 confident genus calls changed, 0 of them explained (a deliberately excluded top hit or a label update), 0 unexplained

No rule fired (see docs/recalibration.md for the rules).

