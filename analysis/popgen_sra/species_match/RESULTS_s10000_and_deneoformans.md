# The small index for the fetch hook, and the two C. deneoformans calls (2026-11-19)

## 1. bfd_s10000 against bfd_s1000 (the hook would load the small one)
Same 2,002 sketches, same code (`popsketch species-match`, defaults), compared with `compare_indexes.py` (`species_match_v3` against `species_match_v3_s10000` on the HPCC).
- Status identical for **1,993 of 2,002 runs (99.6%)**; `exclude_default` identical for **all 2,002**: the 12 runs excluded by default are the same, and every `mixed`, `genus_match`, `host_dominated` and `mismatch` verdict is unchanged.
- The only changes are 9 runs between `match` (5 to `genus_only`) and `genus_only` (4 to `match`), both kept: the small index calls the species slightly less or more often, no panel changes.
- Limit: seven species, one read type; the 9 flips are not examined.

### Cost of one check (SRR30534096, HPCC node, 4 CPUs, 2 repeats each)
| index | size | wall | peak memory |
|---|---|---|---|
| bfd_s10000 | 0.56 GB | 0.6 to 1.5 s | 1.4 GB |
| bfd_s1000 | 5.6 GB | 6.2 to 8.1 s | **13.7 GB** |

With 4 fetches at a time the small index needs about 5.6 GB on top of the fetch processes (about 0.2 GB each); the big one about 55 GB, far over the 8 to 12 GB the P7 and P9 jobs request. The hook uses s10000. Both give `mismatch` for SRR30534096. Not measured: a first read of the index from cold storage on a fresh node (the repeats above are warm).

## 2. SRR34066268 and SRR34066269 (filed C. neoformans, called C. deneoformans by the index)
First 3M read pairs mapped (minimap2 -x sr) to H99 (C. neoformans) and JEC21 (C. deneoformans); the clean C. neoformans controls (DRR829882, DRR829879) map 97 to 98% to H99 and 58 to 60% to JEC21.

| run | strain | to H99 | to JEC21 | reading |
|---|---|---|---|---|
| SRR34066269 | CBS7748 | 58.5% | **96.6%** | the mirror of the controls: *C. deneoformans* |
| SRR34066268 | CBS7816 | 65.1% | 57.9% | no clear parent: like the hybrid-like runs in `../RESULTS_hybrids_and_155.md`, not a clean *C. deneoformans* |

The index's species call for SRR34066268 is not supported by the mapping. `curated_species_proposals.tsv` holds one proposal (SRR34066269) and records why SRR34066268 is not proposed. Nothing was relabelled; the file is for the owner to review. Limits: 3M pairs, one reference per species, a lineage-divergent *C. neoformans* is not excluded for SRR34066268.
