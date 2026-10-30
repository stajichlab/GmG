# Contamination spike-in (#41): bacterial and human DNA in fungal assemblies

Fungal assemblies plus foreign DNA (`scripts/contam_sim.py`, `scripts/hpcc/20_contam.sbatch`): the donor is a share of the total bases, added in
100 kb pieces of *E. coli* + *B. subtilis* (two genomes, ~9 Mb: at most 15% of a 40 Mb assembly without repeating sequence) or of human chromosome 20.
Index 2026.10.2 (with decoys), `known` thresholds. Truth = the host lineage. Assemblies only; contaminated reads are not simulated. Per-genome
results stay on the HPCC (`GmG_runs/contam`, `GmG_runs/contam_species`).

* **Run 1, 24 hosts** (chromosome-level, genus in the reference, species not): 256 inputs, bacteria 1-15%, human 1-30%.
* **Run 2, 44 hosts** whose species is in the reference (post-freeze assemblies, N50 >= 30 kb, 8-70 Mb; 6 genera, 43 of 44 Ascomycota, so this is a narrow
  test): 484 inputs; 32 of them have a correct species call without contamination, and the table below follows those.

## Genus and phylum are unaffected

Run 1: the phylum and genus call of 255 of 256 contaminated inputs equals the call of the host alone, at every fraction (bacteria to 15%, human to 30%).
The one change is a host that was already a divergent, low-confidence Basidiomycota with `distant_non_fungal_hit` before it was contaminated. Run 2:
474 of 474 genus calls unchanged. No `not_fungal` flag in either run: human DNA does not push a fungus onto a decoy.

## ANI falls as (1 - f)^(1/21); species survives to about 10%

Median change of the top-hit ANI against the host alone (run 1): bacteria 1% -0.0006, 5% -0.0026, 10% -0.0052, 15% -0.0076; human 5% -0.0024, 10% -0.0048,
20% -0.0102, 30% -0.0158; the prediction (1 - f)^(1/21) gives 0.0024 / 0.0050 / 0.0107 / 0.0168 for 5 / 10 / 20 / 30%.

**Run 2: species call kept (correct) among the 32 hosts that have a correct species call alone**

| contamination | bacteria | human | median ANI (bacteria / human) |
|---|---|---|---|
| 0 | 32 / 32 | 32 / 32 | 0.9997 |
| 1% | 32 / 32 | 32 / 32 | 0.9993 / 0.9992 |
| 2% | 31 / 32 | - | 0.9987 |
| 5% | 31 / 32 | 31 / 32 | 0.9973 / 0.9975 |
| 10% | 29 / 32 | 30 / 32 | 0.9948 / 0.9952 |
| 15% | 0 / 26 | - | 0.9924 |
| 20% | - | 0 / 32 | 0.9901 |
| 30% | - | 0 / 32 | 0.9845 |

The species cutoff is ANI 0.9938: once the median ANI falls below it the species call is withheld for nearly every input, so there is a cliff, not a slope,
at about 10-15% foreign DNA. **No wrong species call was made at any fraction**: contamination makes the species call missing, never wrong (the genus call
stays). That is the safe direction.

## Species by containment of the reference in the query: not worth it

The open question was whether the species ANI should also use the containment of the reference in the query, which foreign DNA in the query does not change
(sourmash "max containment"). Counted from the columns already in results.tsv, max(ANI of the query, containment of the reference ^ (1/21)) >= 0.9938 would
rescue 10 of the 32 lost species calls at 20% and at 30% human, 4 of 26 at 15% bacteria, 2-3 of 3 at 10%. The reference is itself only partly contained
by a different strain of the species, so its containment is rarely above the cutoff. A third of the calls back, at the price of a rule that is harder to explain, a
risk of a species call for a query that is a superset of a different genome (chimeras, #45), and a second code path in the WASM and CLI. Not adopted. The
limit is documented instead (docs/limitations.md): the species call needs about 90% of the input to be fungal DNA from the species.

## What is not covered

Contaminated reads (a read set is more exposed because the foreign share of reads is often larger), plant and fungal contaminants, contaminants inside a
chimera, and bacterial shares above 15% (not enough distinct bacterial sequence in two genomes). Each is in #41.
