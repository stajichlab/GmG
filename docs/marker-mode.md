# Marker mode: identify a Fusarium isolate from a marker sequence (tef1, rpb2, ...)

GmG's genome indexes keep one k-mer in 1,000-50,000, so a 700 bp Sanger read has less than one hash in them. Marker mode
uses a separate, small index that keeps **every** k-mer of the FusariumID reference sequences, with one reference per
FusariumID strain holding all of that strain's loci. Everything else is ordinary gmg: `gmg run` / `gmg guess`, the same
results table, report and web app.

## Build the index (once; HPCC)

```bash
# 1. parse the FusariumID download (https://www.fusarium.org/page/ID, "sequences in database" zip)
python3 bench/F1_fusariumid/parse_fusariumid.py fusarioidID<DATE>.fas fid/
# 2. one FASTA per strain + gmg metadata
python3 scripts/fusariumid_markers.py prepare fid/records.tsv fid/loci markers/
# 3. sketch at scaled 1 (all k-mers) and build
gmg sketch -k 21 -s 1 -o markers_sk/ markers/refs/*.fa
gmg build --meta markers/meta.tsv --sketches markers_sk/ -o fusariumid_markers.gmgi --taxonomy-version "FusariumID <DATE>"
```

The 2026-06-23 download gives 4,851 strains, 676k distinct k-mers, a 73 MB index; the whole build takes about a minute.
**Do not publish an index built from FusariumID** without the agreement of fusarium.org; it is their curated data.

## Identify

```bash
gmg run reads/ --index fusariumid_markers.gmgi --thresholds markers -o out/
gmg guess --index fusariumid_markers.gmgi --thresholds markers tef1_read.fa
```

A query can be one locus or several (put the loci of one isolate in one FASTA file). `--thresholds markers` is the
calibrated set for this index (species ANI >= 0.995; `data/thresholds/fusariumid_markers.json`). Species of the
*F. oxysporum* species complex are reported as the complex (`species_complex`), as for genomes (`docs/results-schema.md`).
Expect the QC flag `fragmented` (a marker is shorter than 10 kb); it means nothing here.

## How well it works

Held-out test (`bench/F2_marker_mode/`): 20% of FusariumID strains removed from the index and queried locus by locus
(3,408 sequences). See the table there; in short, with `--thresholds markers` about 97% of species calls are right,
species absent from the index are named ~5% of the time, and the share of queries that get a species depends on the locus
(tef1 and rpb2 best, ITS worst). Single loci cannot separate many sister species at all
(`bench/F1_fusariumid/README.md`: half of all species are <= 5 tef1 differences from their closest relative); gmg then
withholds the species (low agreement among the nearest strains) rather than guess. Query several loci for more species calls.
