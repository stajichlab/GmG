# False calls for taxa absent from the reference (BFD benchmark)

`loo` = leave-one-out with relatives present: coverage and precision of the call. `absent` = the genus (or species) was left out of the reference, so any such call is wrong: the false-call rate is the share of those queries that still get a call. Clean labels, 43 suspect genomes dropped.


## genus

| min ANI | leave-one-out coverage | leave-one-out precision | false-call rate when the genus is absent |
|---|---|---|---|
| 0.7765 | 0.986 | 0.969 | 0.579 |
| 0.8000 | 0.975 | 0.974 | 0.433 |
| 0.8500 | 0.941 | 0.981 | 0.265 |
| 0.8800 | 0.921 | 0.983 | 0.208 |
| 0.9000 | 0.906 | 0.985 | 0.183 |
| 0.9200 | 0.889 | 0.986 | 0.132 |
| 0.9500 | 0.862 | 0.989 | 0.096 |
| 0.9700 | 0.830 | 0.990 | 0.084 |

## species

| min ANI | leave-one-out coverage | leave-one-out precision | false-call rate when the species is absent |
|---|---|---|---|
| 0.9700 | 0.841 | 0.973 | 0.247 |
| 0.9800 | 0.821 | 0.975 | 0.208 |
| 0.9900 | 0.773 | 0.978 | 0.152 |
| 0.9938 | 0.722 | 0.982 | 0.117 |
| 0.9970 | 0.626 | 0.985 | 0.068 |
| 0.9990 | 0.446 | 0.987 | 0.024 |
