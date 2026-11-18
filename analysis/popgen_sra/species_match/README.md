# Post-hoc species_match over the 2,002 sketched SRA runs (PopSeq#52; run 2026-11-18 with exact purity and host_dominated)

`popsketch species-match --index bfd_s1000.gmgi -s <set>/sketches --filed filed_<set>.tsv`, filed name = the `species` column of each set's `runs.tsv` (the PopSeq catalog name at the time; the effective-name rule curated > label > NCBI is not applied yet). Defaults `--min-share` 0.02, `--mixed-purity` 0.2, `--host-purity` 0.005. Per-run verdicts are on the HPCC (`GmG_popgen/species_match_v3/sm_<set>.tsv` with a `.meta.json`; `species_match/` is the first run with the top-20-hit purity, `species_match_v2/` the exact-purity run before `host_dominated`); `not_match_or_genus_only.tsv` lists the runs of the v3 run that are neither `match` nor `genus_only`.

Statuses: `match` species agrees; `genus_only` genus agrees, species not called; `genus_match` the index names another species of the genus; `mixed` another genus is called but the filed genus holds at least `--mixed-purity` of the hashes; `mismatch` another genus and little of the filed genus; `host_dominated` no confident call but the filed genus holds `--host-purity` to `--mixed-purity` of the hashes (something else, such as host tissue, dominates the reads: not an isolate genome); `no_fungal_call_suspect` no confident call and almost none of the filed genus (or a large share without a call). **Excluded by default** (column `exclude_default`): `mismatch`, `host_dominated`, `no_fungal_call_suspect`. `purity` is the exact share of the sample's hashes found in any genome of the filed genus.

| set | match | genus_only | genus_match | mixed | mismatch | host_dominated | no_fungal_call_suspect | no sketch |
|---|---|---|---|---|---|---|---|---|
| b_cinerea (300) | 289 | 11 | 0 | 0 | 0 | 0 | 0 | 0 |
| c_beticola (300) | 1 | 294 | 0 | 3 | 0 | 1 | 1 | 0 |
| c_neoformans (300) | 225 | 70 | 2 | 0 | 1 | 1 | 1 | 0 |
| n_discreta (223) | 197 | 26 | 0 | 0 | 0 | 0 | 0 | 0 |
| p_nodorum (300) | 265 | 35 | 0 | 0 | 0 | 0 | 0 | 0 |
| p_oryzae (279) | 258 | 15 | 0 | 0 | 0 | 6 | 0 | 21 (failed fetches) |
| z_tritici (300) | 90 | 209 | 0 | 0 | 0 | 1 | 0 | 0 |

Totals: 1,325 `match`, 660 `genus_only`, 2 `genus_match`, 3 `mixed`, 1 `mismatch`, 9 `host_dominated`, 2 `no_fungal_call_suspect` = 2,002. **12 runs are excluded by default** (1 `mismatch`, 9 `host_dominated`, 2 `no_fungal_call_suspect`).

## What it found
- **SRR30534096 is a `mismatch`** (*Aspergillus*, 27,014 hashes = 33% of the sample; 0.2% of its hashes are in *Cryptococcus*).
- **SRR34066268 and SRR34066269 are `genus_match`**: the index calls *C. deneoformans* (86 to 93% of the hashes in the genus) for runs filed as *C. neoformans*. A species-level mislabel, same study as the runs in `RESULTS_hybrids_and_155.md`.
- **The three *C. beticola* runs (ERR12084585, SRR31166608, SRR31166716) are `mixed`**, not `mismatch`: 32 to 50% of their hashes are *Cercospora* and 47 to 58% sit in *Penicillium* or *Akanthomyces* (5 to 34% of reads map to the *C. beticola* reference). Mixed samples, not mislabels. ERR12084630 (no call, purity 14%) and SRR31166705 (no call, 50%) are suspects of the same kind.
- **`host_dominated` holds nine runs**: the six *P. oryzae* runs of one study (ERR1976666 to ERR1976671, infected wheat, 0.6 to 5.9% of hashes in *Pyricularia*), ERR12084630 (*C. beticola*, 14%), SRR10233430 (5%) and ERR13824414 (*Z. tritici*, 1.4%). The status says the filed genus is present but small, not why: for the wheat study the cause is known; for ERR13824414 (probably non-fungal) and the others it is not. `no_fungal_call_suspect` keeps SRR31166705 (no call, 50% purity) and SRR26462617 (0.1%). All examined in `RESULTS_species_match_suspects.md`.
- `genus_only` (660 runs) mostly measures thin species calls in the index (*C. beticola* 294 of 300, *Z. tritici* 209 of 300), not run problems. Not excluded.

## Calibration of the two cutoffs
- `--min-share` (0.02) only applies to a confident genus call that differs from the filed genus. Every such call here has 33 to 58% of the hashes in the best hit; none was weak. Any value from near 0 to 0.3 gives the same verdicts.
- `--host-purity` (0.005) separates `host_dominated` from `no_fungal_call_suspect` when there is no call. The runs on the two sides are the wheat-tissue runs (0.6 to 5.9%) and ERR13824414 (1.4%) against SRR26462617 (0.1%): there is no gap to put a cutoff in, and the boundary is a placeholder (a run at 1% of the genus is no better understood than one at 0.5%).
- `--mixed-purity` (0.2) separates `mixed` from `mismatch`: filed-genus shares of the other-genus calls are 0.2% (SRR30534096) against 32 to 50% (the three *C. beticola* runs); nothing falls in between. **Any cutoff from about 0.01 to 0.3 gives the same verdicts on these runs**, so the sets do not discriminate it and 0.2 stays as a placeholder. UNMEASURED beyond these seven species and one index (bfd_s1000).

## Known limits
- A `mixed` or `mismatch` verdict says what the reads are, not why (leaf material, a pooled library, a wrong label): that needs the sample metadata.
- The filed name is the set's `runs.tsv` species; the effective-name rule is not applied, and a label change does not yet make a verdict stale.
- One index, one read type, seven species.
