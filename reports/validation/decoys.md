# Non-fungal decoys (#79): what they do and do not do

Index `bfd-2026.10.2` = the BFD reference (23.7k genomes) plus 81 non-fungal eukaryote genomes (`data/decoys/decoys.tsv`: one
genome per class/phylum, round-robin, 43 phyla; animals, plants, algae, oomycetes, amoebae, apicomplexans, parabasalids, ...).
Cost: about 0.4 MB per 100 Mb genome, +1-2% of the index. Evaluated with the binary that has the `not_fungal` ANI gate (>= 0.90).

## Flags

* `not_fungal`: the phylum call or the best hit is a decoy AND the best-hit ANI is >= 0.90 (a real relative among the decoys).
* `distant_non_fungal_hit`: the nearest genomes are decoys but the best-hit ANI is < 0.90. The call is made again without the decoys (added after the truth-set run showed 4.7% of leave-one-out fungi getting a Chordata phylum from a chance decoy hit; #39). No verdict: a divergent fungus and a
  divergent non-fungus both land on whichever decoy is nearest by chance (ANI 0.75-0.88, 3-700 shared hashes).

The gate exists because the first version (no ANI condition) flagged 30 of 1,257 real post-freeze fungi (2.4%) and 10 of 199
IndigoAg assemblies as `not_fungal` (best ANI to a decoy 0.881), and the evidence cannot separate them from true non-fungi at the
same distance (below).

## Results

| set | n | `not_fungal` | `distant_non_fungal_hit` | no flag |
|---|---|---|---|---|
| decoys as queries (each is in the index) | 81 | 81 | 0 | 0 |
| held-out non-fungi that have a decoy relative (Entamoeba histolytica, ANI 0.967) | 1 | 1 | 0 | 0 |
| held-out non-fungi without a close decoy (human, zebrafish, rice, *Toxoplasma*, *Selaginella*, *Thalassiosira*, *P. sojae*, *Trichomonas*) | 8 | 0 | 8 | 0 |
| held-out bacteria and phiX | 3 | 0 | 0 | 3 (`no_match`) |
| post-freeze NCBI fungi (1,257) | 1,257 | **0** | 30 | 1,227 |
| IndigoAg assemblies (199) | 199 | **0** | 10 | 189 |

* The decoys only recognise non-fungi that have a relative among them. Anything farther than ANI 0.90 from every decoy (the
  held-out human/zebrafish/rice at 0.80-0.83 are already "far" from the mouse/maize decoys) gets the soft flag, which says "no
  close relative anywhere", not "not a fungus". Sensitivity for unrelated eukaryotes is therefore low by design; a bigger panel
  helps only for lineages it covers. Bacteria and viruses were never a problem (`no_match`).
* False alarms on real fungi: none. The 30 + 10 distant hits are divergent fungi with no close reference relative; their phylum call
  was already `low_confidence` and equally unreliable before the decoys.
* Why not a lower gate: true non-fungi without a decoy relative (ANI 0.75-0.83, 7-4,000 shared hashes) and divergent real fungi
  (ANI 0.75-0.88, 3-719 shared hashes) overlap on ANI, containment and shared hashes; no threshold separates them.

## Regression on fungal calls (1,257 post-freeze assemblies, old v1 index vs decoy release index)

* status `ok` 1,045 -> 1,045; `low_confidence` 144 -> 147; `no_match` 68 -> 65 (3 inputs moved to a low-confidence call on a decoy).
* genus calls: 1,042 -> 1,042, **0 changed**. Phylum calls changed for 25 inputs, all of which were `low_confidence` before and
  are among the 30 distant hits (the "phylum" is now the nearest decoy's phylum, e.g. Chordata at ANI 0.78: read the flag).
* species calls changed for 8 `ok` inputs (7 *Zygosaccharomyces sapae* -> *Z. rouxii*, 1 *Emmonsia* label with brackets). These are
  taxonomy-normalisation differences between the two release indexes (newer NCBI taxdump), not decoy effects: genus and ANI are identical.

## Consequence for users

`not_fungal` is a reliable positive ("most likely not a fungus"); its absence is not evidence of a fungus. Check
`distant_non_fungal_hit` and `low_confidence` calls with a phylum that is not fungal. See `docs/qc.md` and `docs/limitations.md`.

Per-genome results stay on the HPCC (`GmG_runs/decoys/eval`, `external_big/out_decoy2`).

## Addendum (2026-10-30): the Thalia democratica decoy was removed

The reads validation (#45, reports/validation/reads_contamination.md) found that the *Thalia democratica* decoy (GCF_965202585.1) was the best hit of
*Conidiobolus*-type fungi (Entomophthoromycota) at ANI up to 0.913, closer than any fungal genome in the reference (0.88-0.90): two ZyGoLife assemblies and one
read set got a false `not_fungal` and a Chordata phylum. Hiding that one decoy removed all three false calls. All 81 decoys were also classified with themselves hidden:
only 3 have a best hit at ANI >= 0.85: *Ichthyophthirius* (ciliate) to *Thalia* at 0.909 (not biological: the *Thalia* assembly is the odd one), and *Plasmodium* and
*Dictyostelium* to *Neocallimastix* at 0.90 and 0.89 (AT-rich genomes sharing k-mers; not contamination, so those decoys stay). Release 2026.10.4 has 80 decoys.
