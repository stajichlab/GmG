# Cocci reads 2000000, 2026.10.5

558 inputs scored (0 without a truth row); status: {'ok': 557, 'low_confidence': 1}

Precision = correct among the calls made; called = share of inputs with a call at that rank; recall = correct calls / inputs with truth. Errors are excluded. Truth = the lineage in the truth table (NCBI taxonomy or curated); where it is blank the rank is not scored.

## All inputs

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 558 | 1.000 | 1.000 [0.993, 1.000] | 1.000 |
| class | 558 | 1.000 | 0.996 [0.987, 0.999] | 0.996 |
| order | 558 | 1.000 | 0.996 [0.987, 0.999] | 0.996 |
| family | 558 | 0.998 | 0.995 [0.984, 0.998] | 0.993 |
| genus | 558 | 0.998 | 0.995 [0.984, 0.998] | 0.993 |
| species | 558 | 0.962 | 0.996 [0.987, 0.999] | 0.959 |

## By `n50_bin`


**n50_bin = <10 kb** (n = 558)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 558 | 1.000 | 1.000 [0.993, 1.000] | 1.000 |
| class | 558 | 1.000 | 0.996 [0.987, 0.999] | 0.996 |
| order | 558 | 1.000 | 0.996 [0.987, 0.999] | 0.996 |
| family | 558 | 0.998 | 0.995 [0.984, 0.998] | 0.993 |
| genus | 558 | 0.998 | 0.995 [0.984, 0.998] | 0.993 |
| species | 558 | 0.962 | 0.996 [0.987, 0.999] | 0.959 |

## Most frequent wrong genus calls (3 wrong of 557 called)

| truth | called | n |
|---|---|---|
| Coccidioides | Lecanicillium | 2 |
| Coccidioides | Histoplasma | 1 |

## Most frequent wrong species calls (2 wrong of 537 called)

| truth | called | n |
|---|---|---|
| Coccidioides immitis | Lecanicillium sp. MT-2017a | 1 |
| Coccidioides posadasii | Lecanicillium sp. MT-2017a | 1 |

## Calibration: genus precision by top-hit ANI

| top-hit ANI | n | precision |
|---|---|---|
| [0.95, 0.98) | 9 | 0.889 |
| [0.98, 0.995) | 141 | 1.000 |
| [0.995, 1.01) | 407 | 0.995 |

## Calibration: genus precision by genus vote consistency

| genus vote consistency | n | precision |
|---|---|---|
| [0.95, 1.01) | 557 | 0.995 |

## Calibration: species precision by top-hit ANI

| top-hit ANI | n | precision |
|---|---|---|
| [0.98, 0.995) | 130 | 1.000 |
| [0.995, 1.01) | 407 | 0.995 |

## Calibration: species precision by species vote consistency

| species vote consistency | n | precision |
|---|---|---|
| [0.95, 1.01) | 537 | 0.996 |

**reads, 0.1-1 Gb of sequence** (n = 558)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 558 | 1.000 | 1.000 [0.993, 1.000] | 1.000 |
| class | 558 | 1.000 | 0.996 [0.987, 0.999] | 0.996 |
| order | 558 | 1.000 | 0.996 [0.987, 0.999] | 0.996 |
| family | 558 | 0.998 | 0.995 [0.984, 0.998] | 0.993 |
| genus | 558 | 0.998 | 0.995 [0.984, 0.998] | 0.993 |
| species | 558 | 0.962 | 0.996 [0.987, 0.999] | 0.959 |
