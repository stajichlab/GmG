# False calls for taxa absent from the reference (BFD benchmark)

`loo` = leave-one-out with relatives present: coverage and precision of the call. `absent` = the genus (or species) was left out of the reference, so any such call is wrong: the false-call rate is the share of those queries that still get a call. Clean labels, 72 suspect genomes dropped.


## genus

| min ANI | leave-one-out coverage | leave-one-out precision | false-call rate when the genus is absent |
|---|---|---|---|
| 0.7765 | 0.986 | 0.971 | 0.560 |
| 0.8000 | 0.975 | 0.975 | 0.370 |
| 0.8500 | 0.941 | 0.982 | 0.183 |
| 0.8800 | 0.921 | 0.985 | 0.127 |
| 0.9000 | 0.906 | 0.986 | 0.105 |
| 0.9200 | 0.889 | 0.988 | 0.059 |
| 0.9500 | 0.862 | 0.990 | 0.035 |
| 0.9700 | 0.830 | 0.992 | 0.028 |

## species

| min ANI | leave-one-out coverage | leave-one-out precision | false-call rate when the species is absent |
|---|---|---|---|
| 0.9700 | 0.843 | 0.973 | 0.239 |
| 0.9800 | 0.823 | 0.975 | 0.201 |
| 0.9900 | 0.775 | 0.979 | 0.145 |
| 0.9938 | 0.724 | 0.982 | 0.110 |
| 0.9970 | 0.628 | 0.985 | 0.061 |
| 0.9990 | 0.447 | 0.987 | 0.023 |
