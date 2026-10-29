# Taxonomy: where the lineages come from

Every reference genome carries an **NCBI taxid** (BFD's `NCBI_TAXONID`). `scripts/normalize_taxonomy.py` rebuilds the lineage
from it with a recorded NCBI taxdump (`new_taxdump`, md5 and date go into the index manifest as `taxonomy_version`):

* phylum, class, order, family, genus and species are the NCBI taxa at those ranks, and **the taxid of each rank is stored**,
  so every call in `results.tsv` has a `taxid_<rank>` (blank when NCBI has no taxon at that rank for the lineage, e.g. no
  class for some incertae sedis fungi; GmG skips blank ranks rather than failing them);
* merged taxids are followed; a deleted or unknown taxid keeps the source labels and gets no taxid (`taxid_status`);
* the original source label (the BFD species string) is kept in the `label` column of the metadata table and the index,
  reported as `top_hit_label` next to the NCBI name.

## How different was BFD from NCBI?

For the 23,679 BFD genomes (taxdump of 2026-10-29): 23,671 taxids current, 8 merged, none unresolved. Phylum, class and
order labels were identical; family differs for 6 genomes and genus for 9 (renamed placements, e.g. four *Cladophialophora*
genomes now in *Epibryon*); the species string differs for 224 (0.9%), mostly NCBI's bracketed names for polyphyletic
genera (`[Candida] boidinii`) and quoted placeholder names. So the inconsistency feared in #30 was small; the gain
is stable per-rank taxids and a recorded taxonomy version. Report of the production build:
`taxonomy_normalization.md` in the index release.

## Metadata table

`gmg build --meta` reads a tab-separated table with a header. Required columns: `id taxid phylum class order family genus
species`. Optional: `cluster`, `phylum_taxid class_taxid order_taxid family_taxid genus_taxid species_taxid`, `label`.
Columns are found by name, order does not matter, extra columns are ignored.

```
python3 scripts/normalize_taxonomy.py meta.tsv taxdump/ meta_ncbi.tsv --clusters clusters.tsv --report taxonomy_normalization.md
gmg build --meta meta_ncbi.tsv --sketches sketches/ -o bfd.gmgi --taxonomy-version "NCBI new_taxdump 2026-10-29"
```

Indexes written by this version store the taxids and labels (17 fields per genome); older indexes (9 or 10 fields) still load
and simply have no taxids. Older `gmg` builds cannot read the new indexes.
