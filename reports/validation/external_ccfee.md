# External set: ccfee

Two assemblies of the same strains (AAFTF and shovill): 29 and 38 inputs, 27 strains with both.

## AAFTF assemblies

29 inputs. Status: low_confidence 24, ok 3, no_match 2.

Deepest rank called: family 17, order 5, species 2, (none) 2, phylum 2, genus 1.

| rank | called | share |
|---|---|---|
| phylum | 27 | 0.931 |
| class | 25 | 0.862 |
| order | 25 | 0.862 |
| family | 20 | 0.690 |
| genus | 3 | 0.103 |
| species | 2 | 0.069 |

Phyla called: Ascomycota 26, (none) 2, Basidiomycota 1.

QC flags: short_genome 2, fragmented 2, few_hashes 2, distant_non_fungal_hit 2.

`identical_to_reference` (the input, or a copy, is in the reference): 0 of 29.

## shovill assemblies

38 inputs. Status: low_confidence 33, ok 5.

Deepest rank called: family 23, order 7, genus 3, phylum 3, species 2.

| rank | called | share |
|---|---|---|
| phylum | 38 | 1.000 |
| class | 35 | 0.921 |
| order | 35 | 0.921 |
| family | 28 | 0.737 |
| genus | 5 | 0.132 |
| species | 2 | 0.053 |

Phyla called: Ascomycota 35, Basidiomycota 2, Entomophthoromycota 1.

QC flags: distant_non_fungal_hit 4, long_genome 1.

`identical_to_reference` (the input, or a copy, is in the reference): 0 of 38.

## Agreement of the two assemblies of one strain

| rank | both call | same call | disagree | one only |
|---|---|---|---|---|
| phylum | 27 | 26 | 1 | 0 |
| class | 25 | 24 | 1 | 0 |
| order | 25 | 24 | 1 | 0 |
| family | 20 | 20 | 0 | 0 |
| genus | 3 | 3 | 0 | 0 |
| species | 2 | 2 | 0 | 0 |

Genus disagreements between the assemblers: 0.

Phylum call vs the samples.csv Phylum column, AAFTF: 27 of 27 called agree.

Phylum call vs the samples.csv Phylum column, shovill: 37 of 38 called agree.

