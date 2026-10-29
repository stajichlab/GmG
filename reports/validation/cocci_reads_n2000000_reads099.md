# Cocci reads 2000000 (reads species 0.99)

558 inputs scored (0 without a truth row); status: {'ok': np.int64(557), 'low_confidence': np.int64(1)}

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
