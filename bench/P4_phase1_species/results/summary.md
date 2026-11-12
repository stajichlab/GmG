## Outcome
- 299 result lines for 299 manifest runs: {'ok': 299}; exit codes {'0': 299}
- routes {'aws_odp': 297, 'ena_https': 2}; sampling {'stride': 299}; verify scope {'full': 299}; fallback {False: 299}
- declared concurrency in the sidecars: {20: 64, 24: 235}
- coverage reached: min 16.174x, median 25.000x, 1 runs below 24.9x (runs shorter than the target are kept whole)
- transferred 198 GB for 98 Gb of bases kept from 296 Gb of runs
- `.sra` bytes per base over AWS (297 runs): min 0.36, median 0.42, max 6.78 (the earlier range, n = 3: 0.28 to 0.65)

## Per-fetch cost
| quantity | n | min | median | p90 | max |
|---|---|---|---|---|---|
| fetch total seconds (results.tsv) | n = 299 | min 47.4 | median 448.9 | p90 625.0 | max 4362.2 |
| transfer seconds (sidecar) | n = 299 | min 5.0 | median 21.4 | p90 78.8 | max 2558.9 |
| conversion seconds (sidecar) | n = 299 | min 0.0 | median 50.3 | p90 170.6 | max 854.1 |

Time the sidecar timers do not name (fetch seconds minus transfer, conversion and sketch; per run): median 326 s, 80% of the fetch (p10 55%, p90 85%). It includes the plan lookups, the read counts of the full files, the evenly spaced cut and the checks after it; they were not timed separately, so which of them dominates is UNMEASURED.

Sum of per-fetch seconds: 40.9 task-hours (my estimate from the pilot: 16 to 33)
| fetch seconds per Mb of run | n = 299 | min 0.22 | median 0.74 | p90 0.99 | max 1.89 (estimate used: 0.19 to 0.40)|

## Transfer rate per fetch (MB/s, transfer only, from the sidecars) by how many fetches shared the node
| node | fetches on the node | n | min | median | max |
|---|---|---|---|---|---|
| r12 | 8 | 107 | 1.7 | 18.5 | 30.2 |
| r39 | 8 | 62 | 6.7 | 20.0 | 26.1 |
| r11 | 16 | 128 | 1.9 | 18.9 | 25.6 |

ENA (2 runs): 0.62, 0.81 MB/s

## Jobs (sacct)
| task | node | elapsed | peak RSS GB (of 20 GiB requested) |
|---|---|---|---|
| 31995236_0 | r39 | 01:47:24 | 20.0 GiB |
| 31995236_1 | r39 | 01:55:01 | 10.7 GiB |
| 31995250_2 | r11 | 01:06:35 | 9.3 GiB |
| 31995250_3 | r11 | 01:09:25 | 9.6 GiB |
| 31995250_4 | r11 | 01:07:42 | 9.5 GiB |
| 31995250_5 | r11 | 01:06:59 | 9.1 GiB |
| 31995257_6 | r12 | 00:49:43 | 6.0 GiB |
| 31995257_7 | r12 | 00:49:22 | 5.8 GiB |
| 31995257_8 | r12 | 00:47:58 | 6.7 GiB |
| 31995257_9 | r12 | 00:19:52 | 6.3 GiB |

## Sketches and the sampling check
| sketch seconds | n = 299 | min 10 | median 13 | p90 17 | max 28 |
| sketch files per run (MB) | n = 299 | min 3.1 | median 5.1 | p90 6.2 | max 9.5 |

- `popsketch check`: {'ok': 299}; size ratio 0.50 to 2.35 (median 1.09); coverage ratio 0.04 to 0.92

## The species analysis (popsketch run --collapse mean --collapse-within study)
- audit: {'marked': 125, 'included': 174}; flags {'low_coverage': 107, 'size_high': 16, 'unrelated_to_cohort': 6, 'size_low': 2}
- flagged / total by the three largest studies: PRJEB71961 84/169, PRJNA592329 22/64, PRJEB28420 13/50
- k-mer coverage: low_coverage samples median 2.1 (n = 107), unflagged median 14.0 (n = 174)
- collapse: 299 samples to 122 representatives (largest cluster 68); noise floor in popsketch.json
- label tests: 
  - country: n 66, groups 10, R2 0.055, p 0.8350, strata study
