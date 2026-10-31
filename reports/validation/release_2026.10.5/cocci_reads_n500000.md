# Cocci reads 500000, 2026.10.5

558 inputs scored (0 without a truth row); status: {'ok': 558}

Precision = correct among the calls made; called = share of inputs with a call at that rank; recall = correct calls / inputs with truth. Errors are excluded. Truth = the lineage in the truth table (NCBI taxonomy or curated); where it is blank the rank is not scored.

## All inputs

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 558 | 1.000 | 1.000 [0.993, 1.000] | 1.000 |
| class | 558 | 1.000 | 0.995 [0.984, 0.998] | 0.995 |
| order | 558 | 1.000 | 0.995 [0.984, 0.998] | 0.995 |
| family | 558 | 1.000 | 0.993 [0.982, 0.997] | 0.993 |
| genus | 558 | 1.000 | 0.993 [0.982, 0.997] | 0.993 |
| species | 558 | 0.957 | 0.996 [0.986, 0.999] | 0.953 |

## By `n50_bin`


**n50_bin = <10 kb** (n = 558)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 558 | 1.000 | 1.000 [0.993, 1.000] | 1.000 |
| class | 558 | 1.000 | 0.995 [0.984, 0.998] | 0.995 |
| order | 558 | 1.000 | 0.995 [0.984, 0.998] | 0.995 |
| family | 558 | 1.000 | 0.993 [0.982, 0.997] | 0.993 |
| genus | 558 | 1.000 | 0.993 [0.982, 0.997] | 0.993 |
| species | 558 | 0.957 | 0.996 [0.986, 0.999] | 0.953 |

## Most frequent wrong genus calls (4 wrong of 558 called)

| truth | called | n |
|---|---|---|
| Coccidioides | Lecanicillium | 3 |
| Coccidioides | Histoplasma | 1 |

## Most frequent wrong species calls (2 wrong of 534 called)

| truth | called | n |
|---|---|---|
| Coccidioides immitis | Lecanicillium sp. MT-2017a | 1 |
| Coccidioides posadasii | Lecanicillium sp. MT-2017a | 1 |

## Calibration: genus precision by top-hit ANI

| top-hit ANI | n | precision |
|---|---|---|
| [0.95, 0.98) | 13 | 0.846 |
| [0.98, 0.995) | 76 | 1.000 |
| [0.995, 1.01) | 469 | 0.996 |

## Calibration: genus precision by genus vote consistency

| genus vote consistency | n | precision |
|---|---|---|
| [0.95, 1.01) | 558 | 0.993 |

## Calibration: species precision by top-hit ANI

| top-hit ANI | n | precision |
|---|---|---|
| [0.98, 0.995) | 65 | 1.000 |
| [0.995, 1.01) | 469 | 0.996 |

## Calibration: species precision by species vote consistency

| species vote consistency | n | precision |
|---|---|---|
| [0.95, 1.01) | 534 | 0.996 |

**reads, 0.1-1 Gb of sequence** (n = 107)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 107 | 1.000 | 1.000 [0.965, 1.000] | 1.000 |
| class | 107 | 1.000 | 1.000 [0.965, 1.000] | 1.000 |
| order | 107 | 1.000 | 1.000 [0.965, 1.000] | 1.000 |
| family | 107 | 1.000 | 1.000 [0.965, 1.000] | 1.000 |
| genus | 107 | 1.000 | 1.000 [0.965, 1.000] | 1.000 |
| species | 107 | 0.981 | 1.000 [0.965, 1.000] | 0.981 |

**reads, 10-100 Mb of sequence** (n = 451)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 451 | 1.000 | 1.000 [0.992, 1.000] | 1.000 |
| class | 451 | 1.000 | 0.993 [0.981, 0.998] | 0.993 |
| order | 451 | 1.000 | 0.993 [0.981, 0.998] | 0.993 |
| family | 451 | 1.000 | 0.991 [0.977, 0.997] | 0.991 |
| genus | 451 | 1.000 | 0.991 [0.977, 0.997] | 0.991 |
| species | 451 | 0.951 | 0.995 [0.983, 0.999] | 0.947 |
