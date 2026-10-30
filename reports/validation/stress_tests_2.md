# Second stress test (#41, #45): chimeras, hybrids, organelle/plasmid/virus-only inputs, archaea and more negatives

Index `bfd-2026.10.3-s10000.gmgi` (the reference without the 57 ambiguous-label genomes, with the 81 decoys and the distant-decoy rerun), binary
`a9acf5c` (#93), `known` thresholds, `gmg run --top 10`. Assemblies only. Per-input output stays on the HPCC (`GmG_runs/stress2`); the aggregates, the
scripts (`scripts/stress2_*.py`, `scripts/fetch_small_inputs.py`, `scripts/hpcc/30_stress2.sbatch`) and the pair list
(`data/validation/stress2_pairs.tsv`) are in the repository.

## Summary

* **Chimeras of two same-phylum fungi (72 inputs, 24 genome pairs at 50:50, 70:30, 90:10): every call is the majority parent's taxon, the minority parent's taxon or a rank above;
  none of the 80 inputs got a third taxon.** Phylum and class were correct for 80 of 80 (all inputs including hybrids). At 50:50 the call at the rank where the parents differ is one of
  the two parents (same-class pairs 4 majority + 4 minority; same-order pairs 7 majority, 0 minority, 1 none). At 70:30 and 90:10, of the 32 same-class and same-order inputs 30 were
  called to the majority parent, 1 to the minority parent (a same-class pair at 70:30: order, family and genus all the minority's) and 1 got no call at that rank.
* **The minority is mostly invisible.** The minority parent's taxon is among the top 10 hits for 7 of 8 same-class pairs at 50:50 but for 1-2 of 8 at 30% and 10%; for
  same-order pairs 4/8, 3/8, 3/8. A hit list cannot find a minority below 50%: decomposition (#31) is needed.
* **The top-hit ANI is a mixture meter.** The median ANI falls as p^(1/21) with the majority share p: 0.963-0.968 at 50:50 (0.5^(1/21) = 0.968), 0.977-0.981 at 70:30 (0.983),
  0.989-0.994 at 90:10 (0.995). A top hit of a clearly known genus with ANI 0.96-0.98 and nothing at 0.99+ is what a half-and-half input looks like.
* **Hybrids: simulated (union of two same-genus genomes, 8 pairs) genus correct 8/8, species withheld 8/8**, status `ok`, median ANI 0.969. **Seven real hybrid assemblies with their own species hidden: genus correct 6/7, the seventh kept at the family; with the genus hidden never a genus call (0/7) and the rank above right for 6/7.**
* **Mitochondrion-only inputs (10 RefSeq fungal mtDNAs, 19-135 kb): 7 correct calls (5 to species, 2 to genus), 3 `no_match`, 0 wrong.** They rest on 5-19 hashes and carry
  `few_hashes`. A small-input mode is not worth building now (see below).
* **Plasmid and viruses: no fungal call.** A yeast plasmid (1 hash) and 7 viruses (lambda, T4, mimivirus, PBCV-1, vaccinia, HSV-1, phiX174) all `no_match`.
* **Archaea, bacteria and more non-fungal eukaryotes (13 genomes): 0 confident fungal calls.** 10 `no_match`, 3 `low_confidence` with `distant_non_fungal_hit`.

## Chimeras and hybrids

24 pairs of chromosome-level assemblies (8 to 60 Mb, genus in the reference, same phylum): 8 same genus (different taxid), 8 same order (different genus), 8 same class
(different order). Chimera: majority genome A and minority genome B in 100 kb pieces, amounts cut so the shares are exact; the hybrid is A + B in full. Truth = the
lineage of A; at 50:50 either parent counts as right. "Called (correct)" is at the last rank the parents share; the split columns are at the first rank where they differ.

| kind | variant | % A | n | ok / low / no_match | common rank | called (correct) | at split: A / B / other / none | B in top-10 hits | A in top-10 hits | median ANI |
|---|---|---|---|---|---|---|---|---|---|---|
| same_class | chimera | 90 | 8 | 7 / 1 / 0 | class | 8 (8) | 7 / 0 / 0 / 1 (order) | 1/8 | 8/8 | 0.9933 |
| same_class | chimera | 70 | 8 | 8 / 0 / 0 | class | 8 (8) | 7 / 1 / 0 / 0 (order) | 2/8 | 8/8 | 0.9814 |
| same_class | chimera | 50 | 8 | 3 / 5 / 0 | class | 8 (8) | 4 / 4 / 0 / 0 (order) | 7/8 | 8/8 | 0.9679 |
| same_order | chimera | 90 | 8 | 8 / 0 / 0 | family | 8 (8) | 8 / 0 / 0 / 0 (genus) | 3/8 | 8/8 | 0.9889 |
| same_order | chimera | 70 | 8 | 8 / 0 / 0 | family | 8 (8) | 8 / 0 / 0 / 0 (genus) | 3/8 | 8/8 | 0.9771 |
| same_order | chimera | 50 | 8 | 7 / 1 / 0 | family | 7 (7) | 7 / 0 / 0 / 1 (genus) | 4/8 | 8/8 | 0.9630 |
| same_genus | chimera | 90 | 8 | 7 / 1 / 0 | genus | 7 (7) | species: see below | n/a | 1/8 | 0.9940 |
| same_genus | chimera | 70 | 8 | 8 / 0 / 0 | genus | 8 (8) | species withheld 8/8 | n/a | 1/8 | 0.9811 |
| same_genus | chimera | 50 | 8 | 8 / 0 / 0 | genus | 8 (8) | species withheld 8/8 | n/a | 0/8 | 0.9679 |
| same_genus | hybrid | 100 | 8 | 8 / 0 / 0 | genus | 8 (8) | species withheld 8/8 | n/a | 1/8 | 0.9688 |

Read with care:

* **No wrong call at the common rank in any of the 80 inputs** (`gmg validate` against the majority parent: phylum and class precision 1.000 with every input called). Below the common
  rank `gmg validate` lists order precision 0.937, family 0.973 and genus 0.972; every call it counts as wrong is the *minority parent's* taxon (order 5, family 2, genus 2; none a third
  taxon). 6 of those 9 are 50:50 same-class inputs, where either parent is a right answer; the other 3 are one same-class 70:30 input called to the minority parent. A majority of 70% does
  not guarantee the majority parent is the one called.
* **8 of 80 inputs are `low_confidence`:** 5 same-class 50:50 and 1 same-order 50:50 (no call below the class or order), 1 same-class 90:10 (nothing below the class) and 1 same-genus 90:10
  (family kept, genus not). The 50:50 same-class group is the weakest: 5 of 8 are `low_confidence`.
* **Species:** for same-genus pairs the species is called for 4 of 8 at 90:10 and none at 70:30 and 50:50. The truth species is blank for most of these pairs (new
  assemblies with no species name), so those 4 cannot be scored; none of them was withheld wrongly. The call is withheld whenever the ANI falls under the species cutoff
  (0.9938), which a 70:30 mixture always does (ANI 0.98).
* **Hit lists show the minority only when it is at least half the input**: see the B in top-10 hits column. Below 30% a minority fungus is not distinguishable from
  its absence by `results.tsv` / `hits.tsv`. The reference-side relatives of the majority fill all ten slots.
* **Flags:** no `not_fungal`, no `distant_non_fungal_hit`, no `mixed_lineage` and no other flag than `long_genome` (15 of 80 inputs: the larger assemblies) on any chimera or
  hybrid. The `weak_phylum` flag (assemblies, best-hit ANI < 0.85, #99) is not in the binary used here (`a9acf5c` predates it); counted from the ANI column it would have fired for
  0 of the 80 chimeras and hybrids (the lowest top-hit ANI is 0.894).
* **Pairs:** "same genus" means different taxids of one genus, which includes strains of one species when the species name is missing; the ANI between the parents was not
  measured. A pair of the same species behaves as a plain chimera of one genome.

### What this means for #45 and #31

Two-fungus mixtures at 50:50 to 90:10 never produce a confident wrong taxon: the majority is called, or the rank above it. A mixture is visible only through the ANI (a
known-looking genus with no hit above ~0.97-0.99) and, at 50:50, through a second genus in the top hits; the minority at 30% and below is not recoverable from the present output.
Distinguishing "a divergent strain" from "a mixture" with the same ANI needs the gather-style decomposition of #31 (the share of query hashes explained by the best genome and
what is left). #45 stays open for reads and for mixtures of more than two or of unrelated phyla (the cross-phylum contamination run is in `stress_tests.md` and `contamination.md`).

## Hybrids

The simulated hybrid (both parent genomes in full) is called to the genus of both parents in 8 of 8 pairs, with the species withheld (ANI 0.97: the sketch of a hybrid contains
half of its k-mers from each parent, so the containment ANI against either parent is about 0.97). The status is `ok` and no flag is raised. Real hybrids differ from this union:
their homeologous chromosomes diverge by a few percent, so the k-mers of the two subgenomes partly coincide and the ANI to the nearest parent is higher than for a union.

**Real hybrids** (`scripts/hpcc/31_real_hybrids.sbatch`, index 2026.10.4, same binary; truth = the genus of the parents). Seven NCBI assemblies of known hybrids: *Saccharomyces pastorianus*
(*S. cerevisiae* x *S. eubayanus*), *S. cerevisiae* x *S. kudriavzevii*, *S. eubayanus* (the cold-tolerant parent, not itself a hybrid), *Candida orthopsilosis* (an interspecies hybrid strain; NCBI now
calls the genus *Lodderomyces*), a *Cryptococcus neoformans* x *C. deneoformans* hybrid, *Zygosaccharomyces parabailii* and *Verticillium longisporum*. Six are in the reference
(four are `identical_to_reference`), so the plain run (all 7 `ok`, genus correct for 7/7 under the current names) is a lookup and says little. The test is what happens when the
assembly's own species, and then its own genus, is hidden from the reference (`--exclude-level`, #85):

| assembly | own species hidden: call (best-hit ANI) | own genus hidden: deepest call (best-hit ANI) |
|---|---|---|
| *S. pastorianus* | genus *Saccharomyces* correct; species of another hybrid genome, *S. cerevisiae* x *S. eubayanus* (0.9997) | family Saccharomycetaceae, correct (0.955) |
| *S. cerevisiae* x *S. kudriavzevii* | genus *Saccharomyces* correct (0.9968) | family Saccharomycetaceae, correct (0.976) |
| *S. eubayanus* | genus *Saccharomyces* correct, species withheld (0.9995) | family Saccharomycetaceae, correct (0.916) |
| *Zygosaccharomyces parabailii* | genus *Zygosaccharomyces* correct; species *Z. bailii* (0.9979), a parent species and not the hybrid | family Saccharomycetaceae, correct (0.790, `distant_non_fungal_hit`) |
| *Cryptococcus* hybrid | genus *Cryptococcus* correct; species *C. neoformans* (0.9940), one parent | order Magnaporthales, **wrong** (0.752, `distant_non_fungal_hit`, `low_confidence`) |
| *Verticillium longisporum* | genus *Verticillium* correct; best hit *V. dahliae* (0.9647) | order Glomerellales, correct (0.767, `distant_non_fungal_hit`) |
| *Candida orthopsilosis* | no genus call: family Debaryomycetaceae, best hit *Lodderomyces metapsilosis* (0.838), `low_confidence` | order Serinales, correct (0.796) |

* **Genus correct for 6 of 7 with the species hidden;** the seventh (a close relative at ANI 0.84, below the genus cutoff) is kept at the family, which is right. **Never a genus call with the genus
  hidden** (0 of 7), as it should be, and the rank above is right for 6 of 7 (all `low_confidence`); the one wrong deepest call is the Cryptococcus hybrid, whose nearest remaining relatives
  are at ANI 0.75 (it is called to an unrelated order, `low_confidence` with `distant_non_fungal_hit`).
* **Hybrids are called to the genus of their parents and to a parent species**, never to a third genus. The species call is to one parent's species (or another hybrid genome) at ANI 0.994-0.9997 and
  is withheld otherwise; the sketch cannot say it is a hybrid. This agrees with the simulated union (genus 8/8, species withheld).
* The hybrid assemblies carry no hybrid-specific flag; `high_N` (Verticillium) is the only flag that is not `distant_non_fungal_hit`.

## Mitochondrion-, plasmid- and virus-only inputs

| input | length | status | call | hashes | flags |
|---|---|---|---|---|---|
| *Saccharomyces cerevisiae* NCYC3594 mtDNA | 78,917 | ok | *S. cerevisiae* | 8 | short_genome, few_hashes |
| *Candida albicans* L757 mtDNA | 33,928 | ok | *C. albicans* | 5 | short_genome, few_hashes |
| *Cryptococcus neoformans* var. grubii H99 mtDNA | 24,919 | ok | *C. neoformans* | 5 | short_genome, few_hashes |
| *Agaricus bisporus* H97 mtDNA | 135,005 | ok | *A. bisporus* | 19 | short_genome, few_hashes |
| *Pyricularia oryzae* Guy11 mtDNA | 34,865 | ok | *P. oryzae* | 5 | short_genome, few_hashes |
| *Aspergillus nidulans* FGSC A4 mtDNA | 33,227 | ok | *Aspergillus* | 5 | short_genome, few_hashes |
| *Fusarium graminearum* mtDNA | 95,676 | ok | *Fusarium* | 9 | short_genome, few_hashes |
| *Neurospora crassa* OR74A mtDNA | 64,840 | no_match | | 3 | short_genome, few_hashes |
| *Schizosaccharomyces pombe* MT1 mtDNA | 19,433 | no_match | | 2 | short_genome, few_hashes |
| *Trichoderma reesei* mtDNA | 42,130 | no_match | | 2 | short_genome, few_hashes |
| *Saccharomyces eubayanus* plasmid | 6,161 | no_match | | 1 | short_genome, fragmented, few_hashes |
| phage lambda, T4, phiX174; mimivirus, PBCV-1, vaccinia, HSV-1 (7) | 5 kb to 1.2 Mb | no_match (7/7) | | 1 to 107 | short_genome (few_hashes below ~40 hashes) |

(Rhizopus delemar and Batrachochytrium dendrobatidis mitochondrial genomes were not found by the E-utilities query and are not in the table.)

**Evidence for a small-input mode.** At scaled 10000 a 20-135 kb organelle genome yields 2-19 hashes. The calls that were made rested on 5-19 hashes, were correct 7 of 7
(5 to species, 2 to genus) and were consistent with the `few_hashes` and `short_genome` flags, which are set for every one of them; the 3 inputs with 2-3 hashes returned `no_match`,
which is the right failure (nothing is claimed). No virus, plasmid or bacterium produced a fungal call. Not measured: an organelle genome of a species that is *not* in the
reference (a mitochondrion carries few shared k-mers with a different genus, so the likely result is `no_match`, but it was not tested), and the false-call rate of 5-hash calls in
general; 10 inputs are a sample, not a calibration. **Recommendation: do not build a small-input mode now.** The existing flags tell the user what they have, the failure mode is
`no_match` and not a wrong call, and the only way to get more hashes from a 30 kb input is a finer index (scaled 1000 would make the full index about ten times larger, near 2 GB,
for a niche input class). Revisit if a user needs to identify organelle- or plasmid-only data routinely, and then as a second, optional index of organelles built from RefSeq
mitochondrial genomes (a few hundred MB at most), not as a finer scale of the main index.

## Archaea, bacteria and more non-fungal negatives

Never a confident fungal call: no genus or species call, and none without a flag. The list is `data/validation/negatives_2.txt`.

| genome | status | deepest call | best-hit ANI | flags |
|---|---|---|---|---|
| *Methanocaldococcus jannaschii* (archaeon) | no_match | | 0.830 | short_genome |
| *Sulfolobus acidocaldarius* (archaeon) | no_match | | 0.776 | short_genome |
| *Haloferax volcanii* (archaeon) | no_match | | 0.804 | |
| *Methanosarcina acetivorans* (archaeon) | no_match | | 0.792 | |
| *Cyanidioschyzon merolae* (red alga) | no_match | | 0.744 | |
| *Tetrahymena thermophila* (ciliate) | low_confidence | Neocallimastigaceae | 0.862 | distant_non_fungal_hit |
| *Leishmania major* (kinetoplastid) | low_confidence | Entomophthoraceae | 0.773 | distant_non_fungal_hit |
| *Nematostella vectensis* (cnidarian) | low_confidence | Gigasporaceae | 0.779 | high_N, distant_non_fungal_hit |
| *Naegleria gruberi* (heterolobosean) | no_match | | 0.800 | high_N, distant_non_fungal_hit |
| *Streptomyces coelicolor* | no_match | | 0.776 | distant_non_fungal_hit |
| *Mycoplasma genitalium* | no_match | | 0.851 | short_genome, few_hashes |
| *Pseudomonas aeruginosa* | no_match | | 0.775 | |
| *Wolbachia* (wMel) | no_match | | 0.796 | short_genome |

* The 3 eukaryotes without a close decoy relative get a low-confidence call at a family-level rank from the nearest fungal genomes and `distant_non_fungal_hit`; that flag means
  "no verdict" (docs/qc.md). Two of them (*Leishmania*, *Nematostella*, best-hit ANI 0.77-0.78) are below the 0.85 of the `weak_phylum` rule (#99) and would carry that flag too.
* No `not_fungal` flag appeared on any input of this report, so the decoy-panel issue found elsewhere (a decoy giving a false `not_fungal` on a fungus) does not touch these runs.

## Not covered

Contaminated reads (see the reads reports), real hybrids not in the list above, mixtures of three or more genomes, mixtures across phyla at 50:50, pairs whose parents are
the same species, and an organelle-only input from a species that is not in the reference.
