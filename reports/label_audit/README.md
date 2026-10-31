# BFD label audit

Genomes of the BFD (the lab's own reference, 23.7k fungal genomes) whose label conflicts with their near-identical neighbours, found from sequence similarity only (leave-one-out of the whole BFD, scaled=10000, k=21; `scripts/find_label_suspects.py`, `scripts/analyze_clusters.py`).

| file | what |
|---|---|
| `upstream_report.md` | **Hand-off report (#14)**: every genome or read set whose label looks wrong, by confidence tier (A strong, B moderate, C weak or ambiguous), with a method note, counts and caveats. Written so that a copy can be cut down to accession, current organism, suggested taxon and evidence. Nothing has been sent to anyone. |
| `upstream_report.tsv` | The same rows, machine readable (`scripts/make_upstream_report.py` regenerates both). |
| `exclude_ids.txt` | The 43 genomes of the first audit (phylum/class suspects plus 23 hub references). |
| `reference_exclusions.txt` | The 57 genomes left out of the reference from release 2026.10.3 (the 43 plus 14 found by the truth-set run, #39/#90). |
| `loo_suspects.tsv` | Leave-one-out conflicts at phylum, class or order (ANI >= 0.95, vote consistency >= 0.8). |
| `hub_references.tsv` | 23 genomes that are the nearest genome of >= 3 conflicting phylum/class queries. |
| `size_outliers.tsv` | 69 genomes whose sketch is more than 1.6x the median of their species (kept in the reference). |
| `log_conflicts.tsv` | Class and order conflicts at ANI >= 0.95 from the leave-genus-out run (2,662 rows; mostly genomes whose neighbour is one of the 57). |

All suspicions come from sequence similarity; none has been checked against the material or the submitter.
