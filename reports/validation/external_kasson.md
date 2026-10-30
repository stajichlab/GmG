# External set: kasson

## Classification

36 inputs. Status: ok 16, no_match 16, low_confidence 4.

Deepest rank called: (none) 16, species 12, genus 4, family 3, order 1.

| rank | called | share |
|---|---|---|
| phylum | 20 | 0.556 |
| class | 20 | 0.556 |
| order | 20 | 0.556 |
| family | 19 | 0.528 |
| genus | 16 | 0.444 |
| species | 12 | 0.333 |

Phyla called: (none) 16, Ascomycota 13, Entomophthoromycota 4, Mucoromycota 2, Neocallimastigomycota 1.

QC flags: short_genome 21, fragmented 20, few_hashes 17, identical_to_reference 2, distant_non_fungal_hit 1.

`identical_to_reference` (the input, or a copy, is in the reference): 2 of 36.

## Against the provider names

Truth = the NCBI lineage of the provider genus (names.dmp scientific name or synonym, pinned taxdump); see scripts/score_by_name.py.


status: {'ok': 16, 'no_match': 16, 'low_confidence': 4}; 0 without a usable genus name in the taxdump.

| rank | with truth | called | correct | precision | called share |
|---|---|---|---|---|---|
| phylum | 36 | 20 | 17 | 0.850 | 0.556 |
| class | 36 | 20 | 17 | 0.850 | 0.556 |
| order | 36 | 20 | 16 | 0.800 | 0.556 |
| family | 36 | 19 | 16 | 0.842 | 0.528 |
| genus | 36 | 16 | 14 | 0.875 | 0.444 |

Wrong phylum calls (truth, call): Entomophthoromycota->Mucoromycota x1, Entomophthoromycota->Neocallimastigomycota x1, Ascomycota->Entomophthoromycota x1

Wrong genus calls: 2. Most frequent (provider genus -> call): Cryphonectria->Neonectria, Ophiocordyceps->Massospora


Species, exact name (names with an epithet only): 28 inputs, 12 called, 11 equal to the provider species (0.917 of called).
Species calls that differ: cryphonectria parasitica -> Neonectria punicea.

Metagenome-derived assemblies (herbarium or mixed DNA): 11 inputs, status no_match 5, ok 4, low_confidence 2; genus called for 4.
