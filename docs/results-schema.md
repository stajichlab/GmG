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
| `message` | plain one-sentence reason when `status` is `error` (never a path) |
| `error_code` | `empty_file`, `no_sequence_records`, `not_sequence_file`, `bad_gzip`, `unreadable`, `too_short`, `parse_error`, `internal_error`, `stdin`; blank unless `status` is `error`. `results.jsonl` also has `error_action`, what to do about it |
| `taxid_phylum` .. `taxid_species` | the NCBI taxid of each called taxon (docs/taxonomy.md); blank for a rank that was not called or has no NCBI taxon |
| `top_hit_label` | the original label of the top-hit genome in the source collection, next to the NCBI `top_hit_species` |
| `withheld` | why the first uncalled rank was not called, e.g. `genus: ani_below_cutoff (ANI 0.8870 < 0.9000)`; codes: `no_hits`, `too_few_shared_hashes`, `ani_below_cutoff`, `low_agreement`; blank when a species was called |
| `expected` | the label from the sample sheet, if one was given |
| `expected_match` | `yes` if `expected` equals the call at some rank (case-insensitive), `no` if not, blank without an `expected` |
| `coverage_est` | reads only: estimated genome coverage of the sketched reads (see `qc.coverage_est`), blank otherwise |
| `explained_fraction` | share of the query's hashes found in any reference genome of the genus of the best hit (`qc.mixture.explained_fraction`); blank without a best hit or when the input is `not_fungal` |
| `second_taxon`, `second_explained_fraction` | only when `mixture_suspected` is flagged: the second genus (a different genus from the best hit's) and the share of **all** query hashes that only it explains, i.e. hashes that no genome of the best hit's genus has. About the share of the sample's DNA that is the second organism (assemblies; reads: the share of the k-mers above the abundance cutoff) |
| `schema_version` | 8 |

### Thresholds

* `known` (default, the *balanced* set) was calibrated on BFD for queries whose genus has relatives in the reference (about 95% precision per rank in leave-one-out and leave-species-out tests), with the **genus cutoff raised to ANI 0.90**: with the earlier 0.7765 a genus was named for 58% of queries whose genus is *not* in the reference; at 0.90 it is 18% (`reports/validation/novelty_curve.md`, issue #57; 10.5% on release 2026.10.5 after the 72 genomes with ambiguous labels were left out, `reports/validation/release_2026.10.5/`), for about 8 points of coverage in the common case. Since #58 a **safety margin** applies at genus: the genus vote consistency is 0.9 (was 0.5; the genus ANI stays 0.90), the stricter of the previous cutoff and a fit at 0.97 instead of 0.95; the family cutoff stays ANI 0.7715 (a family margin of 0.8837 was tried and reverted); the other ranks are unchanged (`reports/validation/margin_validation.md`). The thresholds are in the binary, not in the index.
* `permissive` is the earlier set (genus cutoff 0.7765, no tentative genus): more genus calls, many of them wrong when the genus is absent.
* `novel` is far stricter: family, genus and species are never called because no threshold reached 95% precision when the genus was left out.
* Species stay at ANI 0.9938 (a species is still named for 11.7% of queries whose species is absent, 11.0% on release 2026.10.5).
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

`run.json` `thresholds.reads_species_min_ani` (0.99): read sets use this species ANI cutoff instead of the assembly one when the latter is higher (but not when species calls are switched off), because sequencing errors lower the containment ANI of reads.
Read files of one sample are merged into one input and one row: files in the same directory whose names differ only by `_R1`/`_R2` (optionally `_001`), `_1`/`_2`, `.1`/`.2` or a lane token `_L00N`; the row id is the shared stem and `results.jsonl` has `files` (names). `--max-reads` applies to each file. A lone file keeps its own id. `--no-group` turns this off; `--list FILE` adds inputs from a file with one path per line (`#` comments allowed).
`--report-format md` (or `all`) on `gmg run`, and `gmg report DIR --format md`, also write `report.md` and `figures/*.svg` (standalone SVG, colours embedded; renders on GitHub or in a notebook). The Markdown has the same numbers as `report.html`: summary, figures, inputs needing attention, one table row per input (up to 5,000), methods and thresholds. `--report-format pdf` (or `all`) writes `report.pdf` (A4, pure Rust via krilla; figures are vector, text is selectable, Inter is embedded as a subset). Inter is © The Inter Project Authors, SIL Open Font License 1.1 (`crates/gmg-pdf/assets/fonts/Inter-OFL.txt`). The PDF is the Markdown report laid out on pages, so the three formats always agree.
`--sample-sheet FILE`: a tab-separated file with a header and columns `id`, `path` (relative to the sheet) and optionally `expected` and `ref` (the query's own reference genome id, for `--exclude-level`); rows with the same id are the files of one sample (R1/R2, lanes), so pairings no naming rule would find can be stated. `--checksums` records the SHA-256 of every input file in `run.json` (`inputs[].sha256`; slow for large inputs). The progress line shows elapsed time and an ETA.
`gmg guess --format json|tsv|table` prints one JSON object (keys `query, status, best_rank, calls, withheld, tentative_genus, top, hits`), a TSV (header, then one row per query: status, best_rank, the six calls, tentative_genus, withheld, top hit, ANI) or one readable line per query. `-q/--quiet` removes progress and summary lines, `-v/--verbose` adds the withheld reason per input; exit codes are in `gmg --help` (0 finished, 1 could not run, 2 usage).
`qc.mixture` (absent without a best hit or for `not_fungal` inputs): `{primary, explained_fraction, suspected, second}`; `second` is the best other genus, even when it did not pass the rule: `{taxon, rank: "genus", phylum, hashes, explained_fraction, fraction_of_rest, ref_coverage, concentration}`
(`fraction_of_rest` = hashes of the second genus over the hashes the primary genus does not explain; `ref_coverage` = how much of its best genome they cover; `concentration` = share of its hashes that sit in that one genome). `suspected` is the flag `mixture_suspected` (see `docs/qc.md`).
`qc.phylum_split`: up to three `{phylum, hashes}` entries, the query hashes found only in genomes of that phylum (informational; it is too noisy to flag mixtures, see `reports/validation/stress_tests.md`);
`qc.abundance_hist` (reads only): `{counts: [seen 1..12 times], over12, min_abund, median_kept}`, the hashes of a read set by how often they were seen;
`median_kept` is the median abundance of the hashes kept (a crude coverage figure);
`qc.coverage_est` (reads only; number or absent) is the estimated **genome coverage of the reads that were sketched** (so it follows `--max-reads`), and `qc.coverage_basis` says how it was found:
`peak` (the histogram has an error tail, a valley and a peak: the coverage is the mode of the peak) or `poisson` (no peak, coverage below about 3x: the ratio of the numbers of hashes seen 3 and 2 times, else 2 and 1 times, of a Poisson distribution).
The k-mer coverage is converted with the mean read length `L`: `coverage = lambda * L / (L - k + 1)`. No estimate when fewer than 30 hashes are kept or too few were seen twice. For RNA-seq it is only the dominant k-mer abundance, and a heterozygous genome can show a second peak at twice the coverage (the lowest strong peak is reported). Accuracy on simulated and real reads: `reports/validation/reads_calibration.md`.
`--min-abund auto` (reads; also `gmg guess`): the abundance cutoff is taken from the valley of the histogram (the position of the valley, at least 2; 2 when there is no peak) instead of the fixed 2; `qc.abundance_hist.min_abund` records the cutoff used. See `reports/validation/reads_calibration.md` for when it helps and which cutoff is the default.

### Leaving the query out of the reference (validation)

Queries that are themselves in the reference would be found at 100% ANI. To score a call honestly:

- `--exclude-ids FILE` hides the listed reference genome ids (first column, one per line) from every query.
- `--exclude-level self|species|genus` hides, per query, its own genome (sample sheet column `ref`, default the sample `id`), or every reference genome of its species or genus (taken from `ref`'s lineage when it is in the index, otherwise from the `expected` column).

Each affected record in `results.jsonl` gets `excluded: {level, n_genomes}`, and `run.json` records `exclusion` (ids file, number of ids, level). `results.tsv` is unchanged.

### Input from a pipe and subsampling

`gmg guess -` reads one FASTA/FASTQ from standard input (gzip is detected from the bytes; the sample id is `stdin`), e.g. `curl -s URL | gmg guess --index ix.gmgi -`. The global `--max-bases N` stops reading each input file after the record that reaches N bases (a cheap subsample of the start of a file; `--max-reads` counts records). `gmg run` takes `-` too, next to files and directories (once; the row id is `stdin`): `curl -s URL | gmg run - --index ix.gmgi -o out/`. Standard input has no file name, so FASTQ needs `--kind reads` (otherwise it is refused with error code `stdin`); `--checksums` skips it.
