# Contaminated reads (#45, reads half): bacterial and human read pairs in fungal read sets

Real fungal read pairs plus simulated foreign read pairs (`scripts/reads_spike.py`, `scripts/reads_spike_plan.py`, `scripts/hpcc/26_reads_spike.sbatch`).
Hosts (18): 6 *Coccidioides* read sets with a species truth (*C. immitis*, *C. posadasii*; the species is in the reference) and 12 ZyGoLife read sets
whose genus is in the reference (genus truth = the lab label resolved through the NCBI taxdump). Each host contributes its first 300,000 read pairs
(about 2x coverage of a 40 Mb genome) plus M = 300,000 f / (100 - f) pairs simulated from *E. coli* + *B. subtilis* (9 Mb) or human chromosome 20 (64 Mb):
2 x 150 bp, 300 bp fragments, error-free, uniform, both strands. f = 1, 5, 10, 20, 30, 50% of all read pairs. 234 inputs per abundance setting; index
2026.10.3, `known` thresholds, reads species cutoff ANI 0.99. Per-input output stays on the HPCC (`GmG_runs/reads/spike`).

"Kept" means the call is still correct among the hosts whose call is correct without contamination (genus: 14 of 18 hosts; species: 6 of 6 *Coccidioides*).

## Default read settings (minimum k-mer abundance 2)

| donor | foreign pairs | genus kept | species kept | wrong genus calls | median ANI | `not_fungal` | `distant_non_fungal_hit` |
|---|---|---|---|---|---|---|---|
| none | 0% | 14/14 | 6/6 | 0 | 0.9948 | 0 | 0 |
| bacteria | 1% | 14/14 | 6/6 | 0 | 0.9947 | 0 | 0 |
| bacteria | 5% | 14/14 | 6/6 | 0 | 0.9925 | 0 | 0 |
| bacteria | 10% | 14/14 | 2/6 | 0 | 0.9862 | 0 | 0 |
| bacteria | 20% | 14/14 | 0/6 | 0 | 0.9758 | 0 | 0 |
| bacteria | 30% | 12/14 | 0/6 | 0 | 0.9708 | 0 | 0 |
| bacteria | 50% | 12/14 | 0/6 | 0 | 0.9685 | 0 | 0 |
| human | 1% | 14/14 | 6/6 | 0 | 0.9948 | 0 | 0 |
| human | 5% | 14/14 | 6/6 | 0 | 0.9943 | 0 | 0 |
| human | 10% | 14/14 | 6/6 | 0 | 0.9927 | 0 | 0 |
| human | 20% | 14/14 | 2/6 | 0 | 0.9861 | 0 | 0 |
| human | 30% | 13/14 | 0/6 | 0 | 0.9749 | 0 | 0 |
| human | 50% | 11/14 | 0/6 | 0 | 0.9467 | 0 | 0 |

## Minimum abundance 1 (every k-mer counted)

Without contamination only 1 of the 6 *Coccidioides* keeps a correct species call (sequencing errors add k-mers that are in no genome and lower the ANI: median
0.979 against 0.995), so the species rows cannot be followed and the setting is not usable for reads. Genus: 14/14 kept up to 10% bacteria and 5% human, 12/14 from
20% bacteria on and at 10-20% human, 11/14 at 30-50% human; `distant_non_fungal_hit` appears for 2 inputs at 50% human; no wrong genus call. Full table:
`scripts/reads_spike_table.py` on `results_abund1.tsv`.

## What it shows

* **Genus survives to 20% foreign reads and no wrong genus call was made at any share** (0 of 216 contaminated inputs: 18 hosts x 12 variants). The loss at 30-50% is a missing call, not a
  wrong one; the genus is withheld for 1-3 of 14 hosts.
* **The species call is lost earlier than the genus, as for assemblies** (reports/validation/contamination.md): kept at 5% bacterial and 10% human reads, mostly lost
  at 10% bacterial and 20% human, gone at 20% and 30%. ANI falls by about (1 - f)^(1/21) of the foreign share.
* **Human reads hurt less than bacterial reads at the same share, and the reason is the abundance filter.** With minimum abundance 2 a k-mer is kept only if it is seen twice. At 300,000
  foreign pairs the 64 Mb human chromosome is covered at 1.4x at most (10% of pairs: 0.15x), so most of its k-mers are singletons and are dropped. The 9 Mb bacterial
  pair is covered at 1.1x at 10% and 10x at 50%, so its k-mers survive and dilute the containment. The filter therefore removes foreign DNA that is at low coverage
  (a large genome at a small share), which is the usual case for a human contamination of a fungal run, and does nothing for a small genome at high coverage
  (bacteria, plasmids, phage). The effect depends on foreign coverage, not on the share of reads alone.
* **No `not_fungal` flag**: neither donor pushed a fungal read set onto a decoy at the default setting.
* **Not simulated:** sequencing errors in the foreign reads (error-free), adapters, a plant or fungal contaminant, and read sets with a lower host depth than 2x.

## Decision

No change to the thresholds or the reads setting. Documented in docs/limitations.md: a read set tolerates about 10-20% foreign reads for genus and about 5-10% for species; the
default minimum abundance of 2 also hides low-coverage foreign DNA; `--min-abund 1` should not be used for reads.
