# `gmg run` outputs (schema version 2)

```
gmg run GENOMES_DIR [MORE_INPUTS...] --index bfd.gmgi -o results/
```

Inputs are files or directories (searched recursively): FASTA (`.fa .fna .fasta .fas .ffn`) or FASTQ (`.fq .fastq`), optionally gzipped.
The id of an input is its file name without those extensions (two inputs with the same id are an error).
`--kind auto|assembly|reads` (auto = by extension), `--thresholds known|permissive|novel|FILE.json`, `--resume`, `--top N`.

| file | content |
|---|---|
| `results.tsv` | one row per input |
| `hits.tsv` | the top-N reference genomes per input (long format) |
| `run.json` | manifest: tool version and commit, command, times, index (k, scaled, seed, genomes), thresholds, status counts, input list |
| `results.jsonl` | full per-input records (what the TSVs and the report are built from; used by `--resume`) |
| `report.html` | the report: one standalone page (inline SVG, no network), see below. `--no-report` skips it; `gmg report DIR` rebuilds it |
| `relationships.json` | pairwise ANI between the assemblies (2 to 60 of them), for the report heatmap |

Empty cell = not applicable / not called. TSV cells never contain tabs or newlines.

## results.tsv

| column | meaning |
|---|---|
| `query_id`, `path`, `kind` | input id, path, `assembly` or `reads` |
| `status` | `ok` = a call at genus or species; `low_confidence` = calls only at family or above; `no_match` = no rank passed (nothing similar enough in the reference, e.g. not a fungus); `error` = the input could not be processed (see `message`) |
| `best_rank` | deepest rank with a call (`phylum` ... `species`) |
| `call_<rank>` | the taxon called at that rank; empty if withheld. Calls are hierarchical: a rank is withheld if a higher rank failed its threshold |
| `tentative_genus` | a genus the evidence points to but does not support as a call: the genus was not called because ANI is below the genus cutoff (0.90 by default), but ANI is at least 0.7765 and the vote and the higher ranks pass. In that range a named genus is right about 80% of the time even when the genus is in the reference, and wrong when it is not, so it is never a call; the row has status `low_confidence` |
| `conf_<rank>` | vote consistency 0-1: share of the near-best reference hits that agree on the taxon (an evidence value, not a probability) |
| `top_hit`, `top_hit_species`, `top_hit_taxid` | best matching reference genome, its species label and NCBI taxid |
| `ani` | containment ANI of the query in the top hit, `containment_query^(1/k)` |
| `containment_query`, `containment_ref` | fraction of the query's hashes found in the top hit / of the top hit's hashes found in the query |
| `shared` | number of shared sketch hashes with the top hit |
| `cluster` | ANI cluster of the top hit (empty if the index has no cluster column) |
| `n_hashes`, `total_bases`, `n50` | sketch size after filtering, input bases, N50 (assemblies) |
| `qc_flags` | comma-separated, see below |
| `message` | error text |
| `schema_version` | 2 |

### Thresholds

* `known` (default, the *balanced* set) was calibrated on BFD for queries whose genus has relatives in the reference (about 95% precision per rank in leave-one-out and leave-species-out tests), with the **genus cutoff raised to ANI 0.90**: with the earlier 0.7765 a genus was named for 58% of queries whose genus is *not* in the reference; at 0.90 it is 18% (`reports/validation/novelty_curve.md`, issue #57), for about 8 points of coverage in the common case.
* `permissive` is the earlier set (genus cutoff 0.7765, no tentative genus): more genus calls, many of them wrong when the genus is absent.
* `novel` is far stricter: family, genus and species are never called because no threshold reached 95% precision when the genus was left out.
* Species stay at ANI 0.9938 (a species is still named for 11.7% of queries whose species is absent).
Held-out checks: `reports/validation/`. Use `--thresholds FILE.json` for your own (`tentative_genus_min_ani` optional; above 1 switches the tentative field off).

### QC flags (provisional, see issue #29)

`short_genome` (< 5 Mb), `long_genome` (> 150 Mb), `fragmented` (N50 < 10 kb), `high_N` (> 5% non-ACGT), `few_hashes` (< 100 sketch hashes),
`identical_to_reference` (ANI >= 0.9999 and the reference is fully contained). Flags explain a result; they never change a call.

## hits.tsv

`query_id, hit_rank, hit_id, taxid, shared, phylum, class, order, family, genus, species, cluster, schema_version` (the lineage is the metadata of the reference genome).

## report.html

One self-contained page (open it in any browser; no network needed; print it for a PDF until the PDF output exists, issue #23). Sections: summary
(headline percentage and counts), how deep each call went, composition by phylum and genus (colours fixed per phylum), how close the nearest
references are, relationships among inputs (heatmap, only for runs with 2 to 60 assemblies), a sortable/filterable table of all inputs, a detail card per
input (inputs needing attention first: errors, no match, low confidence, then lowest ANI; at most `--report-cards` cards, default 400), and a short
"how to read this" with the exact thresholds. Every chart has a "Table view" with the numbers behind it. Long runs: the table lists the first 5,000 inputs;
all inputs are always in `results.tsv`. The same code renders the report in the browser (WASM `render_report`).

## results.jsonl additions

`qc.abundance_hist` (reads only): `{counts: [seen 1..12 times], over12, min_abund, median_kept}`, the hashes of a read set by how often they were seen;
`median_kept` estimates k-mer coverage.
