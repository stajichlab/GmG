# When to re-calibrate after a reference or taxonomy update (#48)

GmG's cutoffs (ANI per rank, vote consistency, the `not_fungal` gate of 0.90, the species cutoff of 0.9938) were calibrated on one reference.
A new release (more genomes, a newer NCBI taxonomy, genomes removed or decoys added) can move them. The rule is: **build the release, run the
release-stability check, and re-calibrate only if a rule below fires.** The check is `scripts/release_compare.py` (exit code 3 = a rule fired); the
runs are `scripts/hpcc/25_release_stability.sbatch` and `26_release_stability_small.sbatch` (the same binary and `known` thresholds for the old
and the new release on the 1,257 post-freeze assemblies, the IndigoAg assemblies and the held-out non-fungal genomes), the golden panel is
`scripts/golden_check.py`.

## Rules (any one fires the trigger)

| # | rule | why this number |
|---|---|---|
| 1 | **Precision drops.** On the post-freeze set (truth = NCBI lineage), the precision of the calls made at any rank (all calls made, flagged or not: `weak_phylum` calls are kept, so they count) is below the lower end of the 95% Wilson interval of the previous release | a drop outside the noise of the previous measurement; 1,257 inputs give intervals of about +-0.01 at phylum to genus |
| 2 | **Called share moves.** The share of inputs with a genus or species call moves by more than 5 points | silent loss (or gain) of calls changes what users see even if the precision holds. The same applies to the **share of phylum calls flagged `weak_phylum`** (best-hit ANI < 0.85, #99): 5.9% of the post-freeze phylum calls in all four releases, with phylum precision 0.996 without the flag and 0.886 with it; a move of more than 5 points means the reference around the queries changed |
| 3 | **Confident calls churn.** More than 2% of inputs that are `ok` in both releases change their genus call **for a reason other than a deliberately removed top hit or a label update** (`release_compare.py` separates them) | taxonomy renames move a few percent of species names, genus calls should be steadier; more unexplained change means the genome set or the cutoffs moved. A removal changes the calls of the inputs whose nearest genome it was, by design (2026.10.3: 3 of 137 confident IndigoAg genus calls) |
| 4 | **The golden panel changes** (`scripts/golden_check.py` fails) | any change of a call, flag or ANI (> 0.0005) of the 17 pinned genomes must be understood before release |
| 5 | **The reference changes size or composition by more than 10%** (genomes, or clusters at 95%), or the taxonomy dump renames more than 5% of the species labels | cutoffs depend on how dense the reference is around a query |

A fired rule does not mean the thresholds are wrong, it means the check that they are right must be run again: `gmg bench` leave-out calibration on
the new index (`reports/validation/clade_calibration.md` is the previous one), the leave-one/species/genus-out truth-set run (`reports/validation/truth_set.md`)
and the stress tests. Record the decision, with the numbers, in `LOG.md` and in the release notes.

## The genus and family safety margin (#58)

The genus and family cutoffs carry a margin: each is the stricter of its previous value and a fit at precision 0.97 instead of 0.95
(`scripts/margin_calibration.py`, `scripts/hpcc/35_margin_calibration.sbatch`; derivation `reports/validation/margin_derivation.md`; before/after
`reports/validation/margin_validation.md`). Other ranks keep the plain 0.95 fit, species because a margin collapses its coverage (0.55 to 0.11).
When a re-calibration is run, re-derive both: run `margin_calibration.py` with the new bench and the shipped `known` file, and copy the resulting
`min_ani` / `min_consistency` for family and genus only. The cutoffs live in `Thresholds::default()` (crates/gmg-core/src/call.rs) and
`data/thresholds/bfd_s10000.{known,balanced}.json` (a unit test keeps the two equal), not in the index: changing them needs a new binary, not an index rebuild or a
new release asset.

## What does not need re-calibration

* Adding or removing a handful of genomes (the 57 ambiguous-label genomes of 2026.10.3): run the check, expect the churn to be confined to the inputs
  whose nearest genome was removed (listed by `release_compare.py` under "old top hit excluded").
* A taxonomy-only refresh: genus and species labels move, ANI and calls do not; the check separates "same top hit, new label" from real changes.
* Decoys: they add `not_fungal` / `distant_non_fungal_hit` flags; fungal calls must not change (rule 3 catches it if they do).

## Cluster ids

`cluster` is the id of the 95% ANI centroid. Clusters are computed once per BFD build (`gmg cluster`) and carried unchanged into the releases, so the ids
are stable (identical across 2026.10, .1 and .2). Removing genomes can remove a centroid: `scripts/reseed_clusters.py` gives the cluster a new centroid
(the first remaining member in the BFD priority order) so that `cluster` never points at a genome that is not in the index. Re-run `gmg cluster` (and
treat every id as new) only when genomes are added.
