# bfd_meta.tsv refresh, 2026-11-18

Regenerated with `scripts/make_bfd_meta.py` from `Fungi_BFD_runs/tables/species.parquet`
(manifest: 24,363 rows, merge_run_id `taxrebuild-20261118`, built 2026-11-19T05:55:15Z,
sha256 `5845014658a898549f6478c56054238cc047c3c1e249a38ac7664493f6f533af`). Same 8 columns and sort order as before.
Per-field changes: `bfd_meta_changes_20261118.tsv` (old and new value for every changed field).

| | rows |
|---|---|
| old | 23,679 |
| new | 24,363 |
| added ASMIDs | 718 (644 Ascomycota, 61 Basidiomycota, 13 other) |
| dropped ASMIDs | 34 (29 Ascomycota, 5 Basidiomycota) |
| common ASMIDs with a changed field | 33 |

Phylum/class changes (2): GCA_054471075.1 (YJ-2025a) and GCA_054471015.1 (YJ-2025b) are now Basidiomycota
(Cystobasidiomycetes, Microbotryomycetes). The species string still says "Saccharomycotina sp." (the submitted name).
Other changes are genus/family/species renames (Cladophialophora to Epibryon, Lecanicillium to Corniculantispora,
Xylaria to Pseudoxylaria, Mycena to Glaucomycena, blank to Domatiomyces; Cryptococcus gattii VGV to hyracis;
Zygosaccharomyces sapae to rouxii for PYCC10392) and 12 taxid updates.

Where this touches GmG files (not regenerated here): `reports/label_audit/{exclude_ids,reference_exclusions}.txt`
and audit item C014 (YJ-2025b, now correctly Basidiomycota; the exclusion is no longer needed for the phylum reason);
C133/C134 (Buchwaldoboletus B06/B07) are among the 34 dropped ASMIDs; `data/validation/sra_runs.tsv` rows for
SRR24089889 and SRR10915186 carry the old genus/family names.

NOT done: the `.gmgi` indexes. The 718 new ASMIDs have no genome file in `input_clean_genomes` (0 of 718), and the
shared sketch directory is gone, so a rebuild needs a full re-sketch.
