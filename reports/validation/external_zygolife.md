# External set: zygolife

## Classification

552 inputs. Status: ok 394, low_confidence 151, no_match 7.

Deepest rank called: species 204, genus 190, family 108, order 30, class 13, (none) 7.

| rank | called | share |
|---|---|---|
| phylum | 545 | 0.987 |
| class | 545 | 0.987 |
| order | 532 | 0.964 |
| family | 502 | 0.909 |
| genus | 394 | 0.714 |
| species | 204 | 0.370 |

Phyla called: Mucoromycota 520, Entomophthoromycota 7, (none) 7, Glomeromycota 5, Neocallimastigomycota 4, Mortierellomycota 4, Basidiomycota 3, Chordata 2.

QC flags: fragmented 38, distant_non_fungal_hit 34, short_genome 4, not_fungal 2.

`identical_to_reference` (the input, or a copy, is in the reference): 0 of 552.

## Against the provider names

Truth = the NCBI lineage of the provider genus (names.dmp scientific name or synonym, pinned taxdump); see scripts/score_by_name.py.


status: {'low_confidence': 151, 'ok': 394, 'no_match': 7}; 1 without a usable genus name in the taxdump.

| rank | with truth | called | correct | precision | called share |
|---|---|---|---|---|---|
| phylum | 551 | 544 | 526 | 0.967 | 0.987 |
| class | 551 | 544 | 525 | 0.965 | 0.987 |
| order | 551 | 531 | 516 | 0.972 | 0.964 |
| family | 551 | 501 | 459 | 0.916 | 0.909 |
| genus | 551 | 394 | 355 | 0.901 | 0.715 |

Wrong phylum calls (truth, call): Mucoromycota->Glomeromycota x5, Entomophthoromycota->Neocallimastigomycota x4, Entomophthoromycota->Chordata x2, Mucoromycota->Basidiomycota x2, Entomophthoromycota->Mucoromycota x2, Mortierellomycota->Mucoromycota x2, Ascomycota->Mucoromycota x1

Wrong genus calls: 39. Most frequent (provider genus -> call): Capillidium->Batkoa, Conidiobolus->Thalia, Neoconidiobolus->Thalia, Umbelopsis->Rhizopus, Mortierella->Linnemannia, Rhizopus->Circinella, Rhizopus->Mucor, Syncephalastrum->Phycomyces, Phycomyces->Gilbertella, Mortierella->Linnemannia, Syncephalastrum->Pilaira, Mucor->Gilbertella, Cunninghamella->Actinomucor, Ellisomyces->Mucor, Pilaira->Syncephalastrum


Species, exact name (names with an epithet only): 499 inputs, 177 called, 102 equal to the provider species (0.576 of called).
Species calls that differ: rhizopus arrhizus -> Rhizopus delemar; umbelopsis ovata -> Rhizopus delemar; lichtheimia hyalospora -> Lichtheimia corymbifera; circinella muscae -> Circinella umbellata; cunninghamella polymorpha -> Cunninghamella echinulata; mucor circinelloides -> Mucor griseocyanus; mucor circinelloides -> Mucor lusitanicus; mucor circinelloides -> Mucor lusitanicus; rhizopus delemar -> Mucor lusitanicus; rhizopus nodosus -> Rhizopus delemar; rhizopus stolonifer -> Rhizopus delemar; rhizopus stolonifer -> Rhizopus delemar; rhizopus stolonifer -> Rhizopus delemar; rhizopus arrhizus -> Rhizopus microsporus; rhizopus cohnii -> Rhizopus arrhizus; rhizopus lyococcus -> Rhizopus arrhizus; phycomyces blakesleeanus -> Gilbertella persicaria; rhizopus japonicus -> Rhizopus arrhizus; mucor plumbeus -> Gilbertella persicaria; rhizomucor pusillus -> Phycomyces nitens.
