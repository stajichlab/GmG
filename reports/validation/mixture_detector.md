# Mixture detector (#31, narrow form): `mixture_suspected`, `explained_fraction`, `second_taxon`

Decision of 2026-10-30: #31 is delivered as a **flag-only detector** for v0.2, not a decomposer. It never changes the call, never names bacterial or human contamination,
ignores the decoy genomes, honours `--exclude-*` and costs one extra pass over the query hashes. Schema version 8.
Code: `RefIndex::explained_by_genus` (`crates/gmg-core/src/index.rs`), `classify::mixture_scan` (`classify.rs`). Analysis: `scripts/mixture_eval.py`, `scripts/mix_asm.py`;
cluster jobs `scripts/hpcc/35_mixdet_*.sbatch`. Index `bfd-2026.10.4-s10000`, `--thresholds known`. Per-input output stays on the cluster (`GmG_runs/mixdet/out`).

## What it computes

After the final call (and the distant-decoy re-call) the query hashes are split in one pass over the index postings:

1. **Primary genus** = the genus of the best-hit genome. `explained_fraction` = the share of the query hashes found in **any** reference genome of that genus (genus-wide, so conserved and
   shared sequence is credited to the primary genus, not to a second one).
2. Every hash not found in the primary genus is credited, once per genus, to each other genus that has it; the best one is the **second taxon**. `second_explained_fraction` = its hashes over
   **all** query hashes. `results.jsonl` (`qc.mixture.second`) also keeps `hashes`, `fraction_of_rest` (over the hashes the primary genus does not explain), `ref_coverage`, and
   `concentration` (share of the second genus' hashes that sit in its single best genome).
3. **`mixture_suspected`** when the second genus has fraction >= 0.01, hashes >= 20, concentration >= 0.5 **and** fraction of the rest >= 0.3. `second_taxon` / `second_explained_fraction` are filled in
   the TSV only for flagged inputs (the JSON always keeps the best candidate).

Decoys and every genome hidden by `--exclude-ids` / `--exclude-level` are invisible to the pass (the same `skip` closure as the search; tests
`a_second_genus_in_the_input_is_flagged_without_changing_the_call`, `the_mixture_pass_respects_the_exclusion_masks`). Inputs flagged `not_fungal` get no mixture record. A second organism of the
primary's own genus is not looked for.

## Data

* **Pure assemblies (negatives), cluster:** the 674-genome RefSeq truth set with its own genome, species or genus hidden (`35_mixdet_truth.sbatch`: 651 / 641 / 571 inputs have a best hit);
  72 mixture parents alone (genome hidden). Workstation: 108 post-freeze assemblies (not in the reference), the golden panel (17).
* **Contaminated assemblies (should not give a *fungal* flag), cluster:** the 24 stress hosts alone and with *E. coli* + *B. subtilis* at 1-15% (112) or human chr20 at 1-30% (120).
  Workstation: 24 other hosts, the same donors (238 inputs).
* **Assembly mixtures (positives), cluster:** 36 truth-set pairs (12 same genus, 12 same order and different genus, 12 different phyla), B at 0.1, 0.5, 1, 2, 5, 10, 20, 30, 50% of the bases
  (`scripts/mix_asm.py`), both parents hidden from the reference (their relatives stay). Workstation: also 66 pairs of the golden-panel fungi (parents not in the reference), 462 inputs.
* **Read mixtures, cluster:** the real ZyGoLife read mixtures of #45 (32 pairs x 9 shares at 4x, 25 pairs also at 20x; 517 sets, default read settings); 64 and 50 pure sets (A or B alone) at the two depths.
  Workstation: simulated 150 bp reads, 14 pairs, 184 sets.
* The cluster binary was built from this branch before the fraction-of-rest condition was added to the shipped rule; the raw numbers in `qc.mixture` do not depend on it and the rule
  is applied offline by `mixture_eval.py` for every number below.

## How the thresholds were chosen

Every number is pooled false-flag rate over all cluster negatives (2,305 inputs: truth set x3, contamination, parents alone, pure reads) and detection of minorities >= 5% (assemblies, different-genus pairs, 12+12 pairs):

| fraction / hashes / concentration / fraction of rest | false-flag rate | detection (assemblies, >= 5%) |
|---|---|---|
| 0.01 / 20 / none / none | 39.3% | 77% |
| 0.01 / 20 / 0.5 / none | 3.5% | 53% |
| 0.01 / 20 / none / 0.3 | 0.17% | 32% |
| **0.01 / 20 / 0.5 / 0.3 (chosen)** | **0.17%** | **32%** |
| 0.005 / 10 / 0.5 / 0.3 | 0.22% | 32% |
| 0.01 / 50 / 0.5 / 0.3 | 0.17% | 31% |

* A novel genome shares conserved sequence with many genera, so without a rule four inputs in ten have a "second genus". **Fraction of the rest** removes nearly all of it: when the best hit explains only a small
  part of the input (a novel genus: explained fraction 0.02-0.07) a 2% second genus is a 2% share of what is left, while in a mixture of known organisms the second genus *is* what the first leaves over.
  This is the "meaningful share of the remaining hashes" of the spec. **Concentration** (one genome holds the hashes) is the second guard: it alone takes 39% to 3.5%, and one organism is at 0.8-1.0 at every share.
* The 1% / 20-hash floor is about 200 kb of second-genus sequence at scaled 10000. 0.5% / 10 hashes does not detect more at >= 5% and costs a little specificity; the floor only matters at 1-2%.
* On the workstation sets the order of the effects was the same (pure post-freeze assemblies: 8 of 438 pooled negatives flagged with fraction/hashes/concentration alone, 1 with the rest condition).
  A rest condition of 0.2 instead of 0.3 makes little difference (it detects about 7 points more, flags the same genomes plus two post-freeze assemblies whose second genus is a close relative at 5-7%).

## False-mixture rate on pure inputs (rule as shipped)

| negative set | inputs | flagged | rate |
|---|---|---|---|
| truth set assemblies, own genome hidden | 651 | 0 | 0.0% |
| same, species hidden | 641 | 0 | 0.0% |
| same, genus hidden (novel genus) | 571 | 0 | 0.0% |
| mixture parents alone (assemblies, genome hidden) | 72 | 0 | 0.0% |
| contamination hosts alone | 24 | 0 | 0.0% |
| hosts + bacterial DNA 1-15% | 112 | 0 | 0.0% |
| hosts + human DNA 1-30% | 120 | 0 | 0.0% |
| **all pure and contaminated assemblies, cluster** | **2,191** | **0** | **0.0%** (95% upper bound 0.14%) |
| post-freeze assemblies, workstation | 108 | 1 | 0.9% |
| hosts, workstation (24 other genomes) alone / + bacteria / + human | 24 / 118 / 120 | 1 / 5 / 5 | one host (below) |
| real ZyGoLife read sets, one parent alone, 4x | 64 | 2 | 3.1% |
| same, 20x | 50 | 2 | 4.0% |
| simulated read sets, pure (workstation) | 46 | 0 | 0.0% |

The three flagged cases are probably **real mixtures**, not false alarms: (a) two ZyGoLife *Dimargaris* libraries (`p12`) with *Cokeromyces* (Mucorales) explaining 22-36% of their hashes and 10-50%
of a *Cokeromyces* reference genome (*Dimargaris* is a mycoparasite grown on Mucorales hosts); (b) a post-freeze *Metschnikowia* assembly (26 Mb) in which a *Yamadazyma* genus explains 39% of the
hashes and 98% of one *Yamadazyma* reference genome is contained (a co-assembly). Neither is confirmed; both count as false flags above. All the other pure inputs are clean.
**The aim of < 2% on assemblies and reads: met on assemblies (0 of 2,191 cluster inputs; 1 of 217 distinct workstation assemblies); on real reads 2 of 64 libraries (3.1%) are flagged, both probable co-cultures, so
the false rate on reads without them is 0 of 62 (4x) and 0 of 48 (20x).**

Bacterial and human DNA do not produce a fungal second taxon (0 of 232 cluster inputs and no flag beyond the host's own on the workstation).

## Detection versus minority share

Assemblies, cluster, both parents hidden, minority = share of the bases (flagged / inputs, 12 pairs per class):

| pair class | 0.1% | 0.5% | 1% | 2% | 5% | 10% | 20% | 30% | 50% |
|---|---|---|---|---|---|---|---|---|---|
| same genus | 0/12 | 0/12 | 0/12 | 0/12 | 0/12 | 0/12 | 0/12 | 0/12 | 0/12 |
| same order, different genus | 0/12 | 0/12 | 0/12 | 3/12 | 3/12 | 3/12 | 4/12 | 5/12 | 5/12 |
| different phylum | 0/12 | 0/12 | 1/12 | 2/12 | 3/12 | 3/12 | 4/12 | 4/12 | 4/12 |

Only a minority **that the reference explains** can be found. Restricted to the pairs whose two parents are each explained >= 50% by their genus alone (the genus is well represented in the reference):

| pair class (pairs) | 1% | 2% | 5% | 10% | 20% | 30% | 50% |
|---|---|---|---|---|---|---|---|
| same order (3 of 12) | 0/3 | 3/3 | 3/3 | 3/3 | 3/3 | 3/3 | 3/3 |
| different phylum (4 of 12) | 1/4 | 2/4 | 3/4 | 3/4 | 4/4 | 4/4 | 4/4 |

With the minority's genus in the reference the limit is **about 2% of the bases** (the 1% / 20-hash floor; 0.5% to 1% only for a few genomes). A minority whose genus is not in the reference (8 of 12 different-phylum
and 9 of 12 same-order pairs here: the truth-set parents with genus hidden are often thin clades) is not explained and **cannot be flagged at any share, 50% included**; that is a limit of the method, not of the threshold.
A same-genus second organism is not looked for. Workstation, 66 golden-panel pairs (49 different-phylum, parents not in the reference): 2 / 11 / 74 / 100% of the different-phylum inputs flagged at 0.5 / 1 / 2 / 5% of the bases, 100% at 10-50%.

Real read mixtures (ZyGoLife, default settings), flagged / pairs; "4x" = 600,000 pairs, "20x" = 2.7M pairs (the minority is the B parent: share = 100 - A):

| class (pairs at 4x / 20x) | 0.1% | 0.5% | 1% | 5% | 10% | 20% | 50% |
|---|---|---|---|---|---|---|---|
| same species (8 / 6), 4x | 0/8 | 0/8 | 0/8 | 0/8 | 0/8 | 0/8 | 0/8 |
| same genus (8 / 7), 4x | 1/8* | 1/8* | 1/8* | 1/8* | 1/8* | 1/8* | 1/8* |
| same order (8 / 7), 4x | 0/8 | 0/8 | 0/8 | 1/8 | 2/8 | 4/8 | 6/8 |
| different phylum (8 / 5), 4x | 0/8 | 1/8 | 1/8 | 3/8 | 4/8 | 6/8 | 5/8 |
| same order, 20x | 0/7 | 0/7 | 0/7 | 3/7 | 3/7 | 4/7 | 5/7 |
| different phylum, 20x | 0/5 | 1/5 | 1/5 | 1/5 | 4/5 | 4/5 | 3/5 |

\* one pair whose parent is itself flagged (the *Dimargaris*/*Cokeromyces* libraries): flat across shares. Restricted to pairs whose parents are each >= 50% explained and not themselves flagged
(4x): same order 0/5 up to 1%, 1/5 at 5%, 2/5 at 10%, 4/5 at 20%, 5/5 at 50%; different phylum 0/6 at 0.1%, 1/6 at 0.5-1%, 3/6 at 5%, 4/6 at 10%, 6/6 at 20%, 5/6 at 50%. At 20x total
depth: same order 0/4 to 1%, 2/4 at 5-10%, 3/4 at 20%, 4/4 at 50%; different phylum 0/3 to 5%, 3/3 at 10-50% (3 pairs).

**Reads hide the minority through the abundance filter.** `--min-abund 2` drops a k-mer seen once, so a minority covered below about 1x contributes few hashes. The limit is the minority's **own coverage**
(about 0.4-0.8x, i.e. 10-20% of a 4x sample, 5% of a 20x sample), not its share. Simulated reads (workstation, 15 Mb parents) show the same: flagged at 5% at 4x total for 11 of 14 pairs and, at 20x total, at 2% for 9 of 9 pairs (1%: 6 of 9). With
`--min-abund 1` the minority would be visible lower down, but the pure controls lose calls (reads_mixtures.md), so it is not the default.

## Abundance error of the explained fraction

Assemblies, flagged inputs, cluster, median over pairs: the second-taxon fraction is within 0.1-0.6 points of the true minority share up to 5% and within 0.1-1.2 points at 10%, 1.2-2.3 at 20%, 3.4 at 30% and 4.4-7.6
at 50% (same order: 0.0214 at 2%, 0.0463 at 5%, 0.0996 at 10%, 0.198 at 20%, 0.266 at 30%, 0.456 at 50%). The primary's `explained_fraction` is 0.945 / 0.894 / 0.794 / 0.692 / 0.506 for a true 0.95 / 0.90 / 0.80 / 0.70 / 0.50.
It is a **lower bound on the share**: it counts only the hashes the reference explains, so a partly novel minority is under-estimated (different phylum at 50%: 0.42 for a true 0.50).

For reads it is a share of **distinct k-mer content**, not of reads: once the minority is covered at 2x or more its whole genome content counts and the fraction saturates near the genome-size ratio. At 20x total depth
the second fraction is 0.17-0.29 for a true 5-20% minority (error 3-24 points) and the primary's explained fraction 0.44-0.61 where 0.95-0.80 would be the read share. Read it as "this much of the input's sequence
content is explained by another genus", not as a mixing ratio for reads; for assemblies it is a good estimate of the base share when both genera are in the reference.

## Cost

One more pass over the postings of the query hashes, the same cost as `search`. `gmg run` of 108 post-freeze assemblies against the representatives index (`bfd-reps-2026.10.4-s20000`, 4 threads, Apple M-series):
15.80 s on main vs 15.55 s on this branch, user CPU 56.0 s vs 56.5 s (+0.9%): sketching the input dominates (about 0.14 s per input end to end, unchanged from the 0.13 s of the reps index).

## What it does not do

* It does not name a contaminant that is not in the reference (bacteria, human, plants), and does not flag it: those inputs carry `distant_non_fungal_hit` / `not_fungal` or nothing (bacterial and human DNA at 1-30%: 0 flags).
* It does not decompose: no per-taxon calls, no lineage below the genus for the second taxon, no read share. The call is never changed.
* It misses a minority whose genus is not in the reference, a same-genus second organism, and (reads) a minority below about 0.4-0.8x coverage.
* A taxonomic split in the reference (two genera of one clade with shared sequence) or a real co-culture gives a flag on an input that is "one organism" by its label; that is the working definition of what it finds
  (see the *Dimargaris* libraries).

## Decision record

* **Narrow form** (2026-10-30): a flag plus two numbers, because the minority is invisible to the hit list below 30-50% (stress_tests_2.md, reads_mixtures.md) while the hashes the best genus does not explain
  are cheap to attribute; a decomposer would need a model of the abundance filter and of novel minorities.
* The earlier rejected mixture flag (docs/limitations.md: "6% false positives at ~50% sensitivity") used `qc.phylum_split`, hashes found in one phylum only. This one counts hashes explained by a single other genus and
  requires one genome to hold them.
* **A different genus, not family or phylum:** same-family neighbours do produce a "second genus" on single novel organisms (a *Lentinus* assembly with 5% *Sanghuangporus*, a *Myriosclerotinia* assembly with 7% *Sclerotinia*,
  workstation); the fraction-of-the-rest condition removes them without giving up genera of one family, which is where a real mixture of relatives lives.
* **Thresholds** were chosen on the cluster negatives and then checked against the workstation sets (they agree); the detection limits are properties of the rule: if a lower limit is wanted, relax the 1% / 20-hash
  floor (cost: a few percent more flags on novel genomes), never the concentration or the rest condition.
