# Cocci reads 100000, 2026.10.5

558 inputs scored (0 without a truth row); status: {'ok': 556, 'no_match': 2}

Precision = correct among the calls made; called = share of inputs with a call at that rank; recall = correct calls / inputs with truth. Errors are excluded. Truth = the lineage in the truth table (NCBI taxonomy or curated); where it is blank the rank is not scored.

## All inputs

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 558 | 0.996 | 1.000 [0.993, 1.000] | 0.996 |
| class | 558 | 0.996 | 0.996 [0.987, 0.999] | 0.993 |
| order | 558 | 0.996 | 0.996 [0.987, 0.999] | 0.993 |
| family | 558 | 0.996 | 0.995 [0.984, 0.998] | 0.991 |
| genus | 558 | 0.996 | 0.995 [0.984, 0.998] | 0.991 |
| species | 558 | 0.950 | 0.994 [0.983, 0.998] | 0.944 |

## By `n50_bin`


**n50_bin = <10 kb** (n = 558)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 558 | 0.996 | 1.000 [0.993, 1.000] | 0.996 |
| class | 558 | 0.996 | 0.996 [0.987, 0.999] | 0.993 |
| order | 558 | 0.996 | 0.996 [0.987, 0.999] | 0.993 |
| family | 558 | 0.996 | 0.995 [0.984, 0.998] | 0.991 |
| genus | 558 | 0.996 | 0.995 [0.984, 0.998] | 0.991 |
| species | 558 | 0.950 | 0.994 [0.983, 0.998] | 0.944 |

## Most frequent wrong genus calls (3 wrong of 556 called)

| truth | called | n |
|---|---|---|
| Coccidioides | Lecanicillium | 2 |
| Coccidioides | Histoplasma | 1 |

## Most frequent wrong species calls (3 wrong of 530 called)

| truth | called | n |
|---|---|---|
| Coccidioides immitis | Lecanicillium sp. MT-2017a | 1 |
| Coccidioides posadasii | Coccidioides immitis | 1 |
| Coccidioides posadasii | Lecanicillium sp. MT-2017a | 1 |

## Calibration: genus precision by top-hit ANI

| top-hit ANI | n | precision |
|---|---|---|
| [0.9, 0.95) | 2 | 1.000 |
| [0.95, 0.98) | 13 | 0.923 |
| [0.98, 0.995) | 67 | 1.000 |
| [0.995, 1.01) | 474 | 0.996 |

## Calibration: genus precision by genus vote consistency

| genus vote consistency | n | precision |
|---|---|---|
| [0.95, 1.01) | 556 | 0.995 |

## Calibration: species precision by top-hit ANI

| top-hit ANI | n | precision |
|---|---|---|
| [0.98, 0.995) | 56 | 0.982 |
| [0.995, 1.01) | 474 | 0.996 |

## Calibration: species precision by species vote consistency

| species vote consistency | n | precision |
|---|---|---|
| [0.95, 1.01) | 530 | 0.994 |

**reads, 10-100 Mb of sequence** (n = 420)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 420 | 1.000 | 1.000 [0.991, 1.000] | 1.000 |
| class | 420 | 1.000 | 0.995 [0.983, 0.999] | 0.995 |
| order | 420 | 1.000 | 0.995 [0.983, 0.999] | 0.995 |
| family | 420 | 1.000 | 0.993 [0.979, 0.998] | 0.993 |
| genus | 420 | 1.000 | 0.993 [0.979, 0.998] | 0.993 |
| species | 420 | 0.940 | 0.992 [0.978, 0.997] | 0.933 |

**reads, <10 Mb of sequence** (n = 138)

| rank | n | called | precision [95% CI] | recall |
|---|---|---|---|---|
| phylum | 138 | 0.986 | 1.000 [0.973, 1.000] | 0.986 |
| class | 138 | 0.986 | 1.000 [0.973, 1.000] | 0.986 |
| order | 138 | 0.986 | 1.000 [0.973, 1.000] | 0.986 |
| family | 138 | 0.986 | 1.000 [0.973, 1.000] | 0.986 |
| genus | 138 | 0.986 | 1.000 [0.973, 1.000] | 0.986 |
| species | 138 | 0.978 | 1.000 [0.972, 1.000] | 0.978 |
