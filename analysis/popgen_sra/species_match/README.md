# Post-hoc species_match over the 2,002 sketched SRA runs (PopSeq#52, step 2 of the agreed order; run 2026-11-17)

`popsketch species-match --index bfd_s1000.gmgi -s <set>/sketches --filed filed_<set>.tsv`, filed name = the `species` column of each set's `runs.tsv` (the PopSeq catalog name at the time; the effective-name rule curated > label > NCBI is not applied yet). Default `--min-share` 0.02. Full per-run verdicts are on the HPCC (`GmG_popgen/species_match/sm_<set>.tsv` with a `.meta.json`); `not_match_or_genus_match.tsv` lists the runs that are neither.

| set | match | genus_match | mismatch | no_fungal_call_suspect | filed runs with no sketch |
|---|---|---|---|---|---|
| b_cinerea (300) | 289 | 11 | 0 | 0 | 0 |
| c_beticola (300) | 1 | 294 | 3 | 2 | 0 |
| c_neoformans (300) | 225 | 72 | 1 | 2 | 0 |
| n_discreta (223) | 197 | 26 | 0 | 0 | 0 |
| p_nodorum (300) | 265 | 35 | 0 | 0 | 0 |
| p_oryzae (279) | 258 | 15 | 0 | 6 | 21 (the failed fetches) |
| z_tritici (300) | 90 | 209 | 0 | 1 | 0 |

## What it found
- The known cases come out as expected: SRR30534096 `mismatch` (*Aspergillus*, 27,014 hashes, 33% of the sample's hashes); SRR26462617 and ERR13824414 `no_fungal_call_suspect`.
- **Three *C. beticola* runs (ERR12084585, SRR31166608, SRR31166716) are `mismatch`** against *Penicillium* / *Akanthomyces* (47 to 58% of hashes), and ERR12084630 and SRR31166705 are suspects. Only 5 to 34% of their reads map to the *C. beticola* reference (3M pairs). That says they are not clean *C. beticola*, **not** that they are mislabelled: they may be mixed samples (leaf material). Not resolved.
- Six *P. oryzae* runs from one study (ERR1976666 to ERR1976671) are suspects with 0.6 to 5% of hashes in the best hit, and were already off size (1.5 to 2.0 times the assembly).
- `genus_match` is dominated by thin species calls (*C. beticola* 294 of 300, *Z. tritici* 209 of 300): the index calls the genus but not the species, so this status mostly measures the index, not the runs. It is kept and not excluded, as agreed.

## Calibration of `--min-share`
Every `mismatch` has 33 to 58% of the hashes in the best hit; every `no_fungal_call_suspect` has 7% or less. No run fell between them, and no run had a genus call with a share under 2% (the `weak` branch never fired). **Any cutoff from about 0.07 to 0.3 gives the same verdicts on these 2,002 runs**, so these sets do not discriminate between values; 0.02 is kept as the default and the value is UNMEASURED beyond these seven species. It matters only for a sample with a weak call.

## Known limits (to settle before the fetch hook uses it)
- `purity` is computed from the top 20 hits only, so it is a lower bound: it is 0.0000 for the runs whose best hits are all another genus, which does not mean they hold no filed-genus hashes. A proper value needs the filed genus's share over all genomes of the index (an index query, not the hit list).
- A mixed sample whose majority is another genus is called `mismatch` even if the filed species is also present. By the agreed rules purity does not drive the status; this case (the three *C. beticola* runs may be one) needs a decision once purity is exact.
- The filed name here is the set's `runs.tsv` species; the effective-name rule is not applied, and a label change does not yet make a verdict stale.
- One index (bfd_s1000), one read type; checked only for the seven species above.
