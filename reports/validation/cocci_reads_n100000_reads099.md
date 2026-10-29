# Cocci reads 100000 (reads species 0.99)

558 inputs scored (0 without a truth row); status: {'ok': np.int64(556), 'no_match': np.int64(2)}

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
