# RNA-seq reads, species and genus hidden from the reference (#43)

Extends `reports/validation/rna_reads_*.md` (345 species-labelled BFD RNA-seq sets, `scripts/make_rna_truth.py`; each set is cleaned, normalised reads of up to 5 SRA runs of one
species). Those reports score reads of a species that is itself in the reference, so a species call is a lookup. Here the same 345 sets are classified with the sample's species
(leave-species-out) or genus (leave-genus-out) hidden from the reference (`gmg run --sample-sheet --exclude-level species|genus`, #85; the name comes from the `expected` column;
3 of 345 names matched no genome and hid nothing). First 100,000 and 500,000 reads, index 2026.10.3 (with decoys, 57 ambiguous-label genomes left out),
`known` thresholds, the binary with the distant-decoy rerun (#89). Script: `scripts/hpcc/27_rna_exclusion.sbatch`; scoring with `gmg validate` and `scripts/reads_summary.py`. Per-input results on the
HPCC (`GmG_runs/reads/rna_exclusion`).

Cells: *called share / precision (calls made)*.

## All calls (500,000 reads)

| exclusion | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|
| none (baseline) | 1.00 / 0.991 (344) | 1.00 / 0.985 (336) | 0.99 / 0.988 (332) | 0.99 / 0.991 (330) | 0.93 / 0.997 (319) | 0.68 / 0.962 (236) |
| species hidden | 0.95 / 0.914 (327) | 0.95 / 0.909 (319) | 0.93 / 0.909 (309) | 0.86 / 0.869 (289) | 0.47 / 0.875 (160) | 0.10 / 0.000 (33) |
| genus hidden | 0.90 / 0.842 (310) | 0.89 / 0.800 (300) | 0.82 / 0.770 (274) | 0.68 / 0.626 (227) | 0.12 / 0.000 (40) | 0.03 / 0.182 (11) |

## Confident calls only (status `ok`, 500,000 reads)

| exclusion | n ok | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| none | 321 | 1.00 / 1.000 (321) | 1.00 / 0.997 (317) | 1.00 / 0.997 (315) | 1.00 / 0.997 (313) | 1.00 / 0.997 (319) | 0.74 / 0.962 (236) |
| species hidden | 164 | 1.00 / 0.994 (164) | 1.00 / 0.988 (164) | 0.99 / 0.988 (163) | 0.99 / 0.950 (161) | 0.98 / 0.875 (160) | 0.20 / 0.000 (33) |
| genus hidden | 45 | 1.00 / 0.956 (45) | 1.00 / 0.889 (45) | 0.98 / 0.864 (44) | 0.98 / 0.691 (42) | 0.93 / 0.000 (40) | 0.24 / 0.182 (11) |

At 100,000 reads the numbers are the same in kind: baseline genus 0.997 (called 0.93), species 0.962 (called 0.68); species hidden: genus called 0.47 (0.887), species called 0.10 (precision 0);
genus hidden: genus called 0.11. The called share of the higher ranks is lower than at 500k (phylum 0.74 species hidden, 0.46 genus hidden) because fewer inputs reach a call.

## What it shows

* **With nothing hidden the earlier result holds**: genus 0.997 (called 0.93) and species 0.962 (called 0.68) at 100k and 500k reads; phylum to genus of the `ok` calls 1.000-0.997.
* **Species hidden (a new species of a known genus):** half of the inputs keep a genus call (0.47 called, precision 0.875; the 20 wrong genus calls were not reviewed one by one: a label conflict or a close genus is likely, but that is not checked here), and 33 inputs (9.6%) get a species call that is necessarily wrong. For assemblies the same test gave 6.9% wrong species calls (reports/validation/truth_set.md): reads are
  slightly more prone to it, since a transcript set has fewer informative k-mers and the species cutoff for reads is 0.99, not 0.9938.
* **Genus hidden (a new genus of a known family):** a genus is still called for 40 inputs (11.6%, all wrong by construction; 12% for assemblies), and 11 get a species call.
  Family precision falls to 0.63 over all calls and 0.69 over the confident ones, against 0.62 for assemblies in the same test.
* **The precision of all phylum calls falls to 0.91 / 0.84 (species / genus hidden) because most of the extra calls are `low_confidence`**: of 28 wrong phylum calls with the
  species hidden, 27 are `low_confidence`; of 49 with the genus hidden, 47. The `ok` calls stay at 0.994 (species hidden) and 0.956 (genus hidden), against 0.997 and 0.966 for assemblies.
  Read the status: a `low_confidence` phylum is a guess (0.82-0.83 precision), an `ok` call is not.
* **Decoys:** 33 (species hidden) and 80 (genus hidden) inputs carry a decoy flag; 32 and 79 are `distant_non_fungal_hit` and, after the distant-decoy rerun, none of them gets a
  non-fungal phylum call. The one `not_fungal` (in every run, also with nothing hidden) is *Puccinia sorghi*, a maize rust: its top hit is the maize decoy (*Zea mays*) at ANI 1.000,
  because RNA of an infected leaf is largely host RNA. That flag is correct: the set is mostly maize, and the call is a non-fungal phylum (Streptophyta), which is what the flag exists for.

## Acceptance

The existing reports plus this one cover RNA-seq at 100k, 500k and 2M reads (the 2M baseline is in `rna_reads_n2000000_*.md`) with the species and the genus hidden: genus precision
0.997, species 0.962 (called 0.68) in the lookup case, and the honest figures for new species and genera above. Not covered: single-sample RNA-seq with strand-specific or total
RNA (rRNA-heavy) libraries, and reads from the spliced-transcript-aware reference (#33, not built).
