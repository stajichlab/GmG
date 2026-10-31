# Tutorial: from a folder of genomes (or reads) to a report

You need `gmg` (a release binary, the container `ghcr.io/stajichlab/gmg`, or `pixi run build`) and a reference index
(`bfd-*.gmgi`; see [index-distribution.md](index-distribution.md)). `gmg info INDEX` shows what an index is made of.

No genomes at hand? `scripts/fetch_examples.sh examples_data` downloads six small public assemblies from NCBI (see
[../examples/README.md](../examples/README.md)); use `examples_data/` below.

## 1. Genomes

```
gmg run my_genomes/ --index bfd-2026.10-s10000.gmgi -o results/
```

`my_genomes/` is searched recursively for `.fa .fna .fasta .fas .ffn` (optionally `.gz`, `.zst`, `.bz2` or `.xz`; the browser app reads only `.gz`). While it runs you see one line
per input with the elapsed time and an ETA. Interrupted? Add `--resume` and finished inputs are reused.

`results/` then holds:

| file | what |
|---|---|
| `results.tsv` | one row per input: status, deepest rank called, the call at each rank, confidence, top hit, ANI, QC flags (columns: [results-schema.md](results-schema.md)) |
| `hits.tsv` | the top reference hits of every input (long format) |
| `report.html` | the report: summary figures, the table, a card per input with the evidence (open it in a browser) |
| `run.json` | tool version, index, thresholds, timing: how the run was made |
| `results.jsonl` | the raw records the TSV and report are made from |

Read `status` first: `ok` (a genus or species was called), `low_confidence` (a higher rank only; the `tentative_genus`
column names the nearest genus when the evidence stops short of a call), `no_match` (nothing in the reference is close
enough: not a fungus, a lineage the reference lacks, or too little of a genome) and `error` (the file could not be read;
`error_code` and `message` say why).

```
column -t -s$'\t' results/results.tsv | less -S
```

One genome without a folder: `gmg guess --index INDEX genome.fa.gz --format table` prints a line such as
`genome.fa.gz  Agaricus  (genus, ANI 0.9846; species withheld: ANI 0.9846 < 0.9938)`; `--format tsv` and `--format json`
give the same information for scripts. The `withheld` column of `results.tsv` says the same for every input: the rank that
was not called and the number that stopped it.

## 2. Reads (DNA or RNA)

```
gmg run reads/ --index bfd-2026.10-s10000.gmgi -o results/ --max-reads 1000000
```

FASTQ files are recognised by extension. R1/R2 and lane-split files of one sample (`S1_L001_R1_001.fastq.gz`,
`S1_L001_R2_001.fastq.gz`, ...) are merged into one row. Pairings that names do not show can be given in a sample sheet:

```
id	path	expected
strainA	data/forward.fastq.gz	Fusarium gamma
strainA	data/reverse.fastq.gz
strainB	data/assembly.fa
```

```
gmg run --sample-sheet samples.tsv --index bfd-2026.10-s10000.gmgi -o results/
```

`expected` is optional; with it, `results.tsv` gets `expected_match` (yes/no) so you can see disagreements at a glance.
Reads use a species cutoff of ANI 0.99 (errors lower the ANI of a read set); about 0.5x coverage (100,000 reads of a
30 Mb genome) is enough for a genus call.

**Long reads (ONT, PacBio).** The same command works. HiFi and ONT R10.4 (error <= 3%) need nothing special (add `--min-abund auto` at 10x or more); for ONT R9.4 and PacBio CLR (5-15% error) expect a genus, rarely a species, from about 5x, and
look for the `high_error_reads` flag in `qc_flags`. Never use `--min-abund 1` for long reads. Details and numbers: `reports/validation/long_reads.md`.

## 3. Other report formats

```
gmg report results/ --format all     # report.html, report.md + figures/*.svg, report.pdf
```

or `--report-format all` on `gmg run`. All three come from the same run.

## 4. How far to trust a call

* The default thresholds name a family at ANI >= 0.7715, a genus only at ANI >= 0.90 with a vote consistency of at least 0.9, and a species at ANI >= 0.9938 (0.99 for reads); they were
  calibrated on the BFD reference (`reports/validation/README.md`). A genus that is not in the reference is usually
  reported as `low_confidence` with a `tentative_genus`, not as a call.
* `--thresholds novel` is the stricter set for queries that may be a new genus; `--thresholds permissive` is the older,
  lower genus cutoff.
* QC flags (`short_genome`, `fragmented`, `few_hashes`, ...) explain a result; they never change a call.
* Contaminated or mixed samples are classified as one organism. The flag `mixture_suspected` (with `second_taxon` and the `explained_fraction` columns) marks an input in which a second fungal *genus* that is in the reference explains at least 1% of the hashes; it cannot see a second organism that is not in the reference, bacteria or human DNA, or one of the same genus (`docs/qc.md`).
* GmG assumes the input is a fungus. Bacterial and viral genomes come out `no_match`, but plant, animal and protist genomes can be given a (wrong) fungal phylum, class, order or family at low confidence; they never reach a genus or species call (`reports/validation/negative_controls.md`).

## 5. Shell completions and manual page

```
gmg completions bash > ~/.local/share/bash-completion/completions/gmg
gmg completions zsh  > ~/.zfunc/_gmg
gmg man | gzip > ~/.local/share/man/man1/gmg.1.gz
```

Exit codes: 0 finished (unreadable inputs are rows with `status=error`), 1 could not run, 2 usage error.

## Checking accuracy against known answers

With a table of what the inputs really are, `gmg validate` reports accuracy per rank and per clade:

```
gmg validate results/results.tsv truth.tsv -o validation/ --title "My test set"
```

`truth.tsv` has an `id` column (the `query_id` of `results.tsv`), any of `phylum class order family genus species` (blank = unknown,
not scored) and optional grouping columns. `validation/validation.html` shows precision and the share called per rank (with 95%
intervals), calibration (precision by ANI), the most frequent wrong calls and the same tables by group; `validation.md` and the
long-format `validation.tsv` hold the numbers. Precision is correct among the calls made; recall counts an input with no call as missed.
