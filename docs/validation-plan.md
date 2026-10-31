# Validation plan (issues #39-#46)

What is already measured, what each remaining validation needs, and in which order to do it. Numbers are from
`reports/validation/` and the NCBI metadata of 2026-10-29. "Effort": S = a day or less, M = a few days, L = a week or more.

## What a validation here means

Every validation answers the same question in the same format, so results can be compared and a regression is visible:

* **Truth table** (`data/validation/<dataset>.truth.tsv`): `id, phylum, class, order, family, genus, species` (NCBI names), plus
  `source` (where the label comes from), `evidence` (`type_material`, `reference_genome`, `curated`, `submitter`, `assembly_label`),
  `in_bfd` (`exact`, `species`, `genus`, `none`: what the reference holds for it) and free `notes`.
* **Run** with the production index and the default thresholds, and with `--thresholds novel` where the dataset can contain new taxa.
* **Report** from `gmg validate RESULTS.tsv TRUTH.tsv -o DIR` (#25; writes validation.md, validation.tsv and validation.html with figures; `scripts/validate_run.py` is the reference implementation it was checked against: identical numbers on 132 table rows of 1,257 results): per rank precision with a Wilson 95% interval,
  fraction called, recall, break-downs by clade, assembly level / depth, and whether the taxon is in the reference;
  the confusion list (top 20 wrong calls); run time and memory.
* **Rules:** a call is *wrong* only when it disagrees with the truth at that rank; a withheld rank is not an error (it lowers
  *called*); a dataset whose truth is blank at a rank is not scored there; queries whose taxon is absent from the reference are
  reported separately (they cannot be right, the question is whether the tool refuses: the false-call rate).

## Shared prerequisites (do these first; they unblock four of the eight)

1. **Reference exclusion in `gmg run` (#85, S-M).** Most truth data is *in* BFD (652 of the 674 RefSeq fungal genomes; 257
   of 267 type-material assemblies), so an honest test removes the query from the reference. `gmg bench` can do leave-one-out,
   leave-species-out and leave-genus-out for genomes already sketched, but `gmg run` cannot, and reads never could. Add
   `--exclude-ids FILE` (ids to drop from the index at search time; the code already has a `skip` closure in `search`) and
   `--exclude-level self|species|genus` driven by the truth table, so one run mimics "this genome/species/genus is not in the
   reference". Needed by #39, #42, #43 (strain-matched reads) and #44.
2. **Dataset registry** `data/validation/datasets.tsv`: name, n, source, truth evidence, license, local path on the HPCC, which index it
   was run against, report path. A table, not code; keeps "which sets exist" in one place (acceptance of #40 is partly counting).
3. **`gmg validate` / validation dashboard (#25, M).** The same report engine (HTML/Markdown/PDF), from `results.tsv` + truth table, so
   each dataset below ends in a page rather than a Markdown by hand.
4. **Disagreement log** `data/validation/disagreements.tsv`: `dataset, id, rank, truth, call, category, evidence, decision, date`;
   categories: `tool_error`, `reference_label_error`, `truth_label_error`, `species_complex`, `taxonomy_synonym`, `contamination`,
   `low_quality_input`. Used by #39 first, then by every other set.

## #39 Trusted truth set and manual review

**Why.** BFD labels have errors (`reports/label_audit/`: 43 suspect genomes, 23 hubs, 69 size outliers), so accuracy measured against
BFD labels mixes tool errors with label errors.

**Data.** NCBI RefSeq fungal genomes: **674** (668 reference/representative, **267 from type material**, levels: 68 chromosome,
69 complete, 324 scaffold, 213 contig). 652 are in BFD; 10 type-material assemblies are not (they are a small fully independent
set). Add curated lab sets (Coccidioides, MAT project) as `curated` evidence. Get the list with
`datasets summary genome taxon fungi --assembly-source RefSeq --as-json-lines` (the `type_material` field).

**Protocol.**
1. Truth = NCBI lineage of the type-material / reference accession (taxdump pinned, as for the index; `scripts/normalize_taxonomy.py`).
2. Evaluate three ways with prerequisite 1: *leave-one-out* (the genome removed: "a new strain of a known species"), *leave-species-out*
   (a new species of a known genus), *leave-genus-out* (a new genus of a known family). For the 10 outside BFD, plain runs.
3. Report each separately: precision/called by rank, species-level recall for leave-one-out (does the type strain find its species?),
   the false-call rate for leave-species/genus-out (does it refuse?).
4. **Review every disagreement by hand** (expected a few dozen): look at the top hits' labels, the cluster of the top hit
   (`gmg` already gives the ANI cluster), the label audit lists, the taxonomy history of the names (synonyms are not errors) and
   record the verdict in the disagreement log. Where the truth label is wrong, say so upstream (#14).

**Deliverable.** `reports/validation/truth_set.md`: accuracy on the trusted set next to accuracy on the full BFD benchmark (the
gap is the label-error share), and the disagreement table. **Acceptance:** >= 250 type-material genomes scored in all three
modes; every disagreement categorised; the report states the share of BFD-benchmark errors that were label errors.
**Effort:** M (S once prerequisite 1 exists; the manual review is the time). **Needs from you:** none, except a view on the review
(I propose to do the first pass and list the uncertain ones for you).

## #40 Independent assembled-genome datasets across clades

**Done.** Post-freeze NCBI set: 1,257 genomes, no quality filter, 0 errors; precision phylum .990, class .991, order .986,
family .983, genus .988 (balanced default) with genus called 84%; absent-genus false names 3/22; N50 < 10 kb assemblies mostly
unresolved (`post_freeze_balanced.md`). That is one dataset (>= 500 outside BFD).

**To add (target 6+ datasets, >= 500 genomes outside BFD, per-clade results):**

| dataset | source and size | truth | what it tests | effort |
|---|---|---|---|---|
| D2 NCBI fungal MAGs | `datasets summary genome taxon fungi --mag only`: **985 MAGs** (lichen genera Umbilicaria, Cladonia, Debaryomyces, Geotrichum; most released 2023-2026) | submitter taxonomy (`submitter` evidence): coarse; report to genus | bins of mixed quality and contamination; lichen/symbiont fungi | S |
| D3 ZygoLife / thin clades | NCBI accessions from the ZygoLife projects (Mucoromycota, Zoopagomycota, Kickxellomycota); chytrids and Microsporidia from NCBI | NCBI taxonomy | clades where BFD has few genomes: the clade-out weakness (#58) | S-M |
| D4 JGI Mycocosm | portal downloads (login needed for bulk); a few thousand genomes, many also in NCBI (count to be checked) | JGI taxonomy + NCBI | access first: ask whether an account/export exists, otherwise skip | M, access-dependent |
| D5 lab assemblies | Coccidioides pangenome, MAT-project genomes, the 14 desert isolates (existing mash IDs for comparison) | lab/curated | closest-to-use case; comparison with the earlier mash IDs | S once the files are located |
| D6 pathogen panels | by genus from NCBI, not in BFD: Aspergillus (1,759 assemblies on NCBI), Fusarium (2,250), Candida (321), Cryptococcus (285), Trichophyton (200) | NCBI taxonomy, type-material subset flagged | species-complex handling inside medically relevant genera; the species the tool must not confuse | S |
| D7 type-material outside BFD | the 10 from #39 | type material | tiny but gold-standard | S |

Per dataset: the standard report, per-clade rows, and the confusion list. Cross-dataset table in `reports/validation/README.md`.
**Needs from you:** where the lab assemblies live; whether a JGI export exists.

## #41 Stress tests

**Done (`reports/validation/stress_tests.md`).** Completeness 10-75% (24 genomes): genus calls unchanged down to 10%
(`short_genome` fires at <= 25%); 5 kb fragmentation: no effect on calls, `fragmented` fires; cross-phylum contamination: 5% flips
the phylum for 3/24 (divergent hosts whose own reference is distant), 20% for 5/24, 50% for 7/24; negative controls: random and
repeat sequence `no_match`; detection of mixtures from phylum-exclusive hashes measured and rejected (#64);
real non-fungal genomes (#26/#79).

**Remaining, with protocols:**
* *Bacterial/human contamination* at 1-30% (E. coli, B. subtilis, human are on the HPCC now): same simulator
  (`scripts/stress_sim.py` + a donor choice), expect hashes of the donor to be absent from the fungal reference, so calls should be
  stable; measure the species-ANI loss `(1-f)^(1/21)` and the point where species calls drop (the earlier review's open question:
  should species use containment of the reference in the query?). S.
* *Same-phylum chimeras* (two fungi, 50:50 to 90:10, same genus / same order / different class): which taxon is reported, does
  consistency fall. S. The answer feeds #45.
* *Hybrids and diploids:* simulate (F1 hybrid = both parents' sequence), and real ones (Saccharomyces pastorianus lager yeast,
  S. cerevisiae x S. kudriavzevii wine hybrids, Candida hybrids, Cryptococcus AD hybrids) from NCBI. Expected: calls to the closer
  parent with lowered consistency at species; document. S-M.
* *Mitochondrion-only and plasmid-only inputs* (NCBI RefSeq organelle genomes of fungi; 20-80 kb: a handful of hashes at scaled
  10000): expected `few_hashes`/`no_match`. Decide whether a **small-input mode** (smaller `scaled`, e.g. 100, for inputs under
  1 Mb) is worth building: it needs a second index at fine scaling (the 5.6 GB s=1000 index is too big; a representatives-only
  s=1000 index would be roughly 0.5-1 GB: an estimate to measure). M, decision needed after measuring.
* *Negative controls* (bacteria, plants, human, viruses must be `no_match`): the plants/animals/oomycetes now answer "not a fungus"
  through decoys (#79); add human and a few archaea/viruses to the held-out list. S.

**Deliverable:** one table: gradient x call correctness x flags; "where it degrades" per axis. **Acceptance:** every bullet above
has a row in `stress_tests.md`.

## #42 DNA reads

**Done.** 558 Coccidioides R1 read sets at 100k/500k/2M reads: genus called >= 99.6%, correct 99.1-99.3%; species (immitis vs
posadasii) called 95-96% at precision .994-.996 after the reads species cutoff (#32), including a leave-immitis-out check (0.6% false
species).

**To add.**
1. **Other genera, other labs.** Harvest runs with trustworthy labels: for BFD genomes with a BioSample, query ENA
   (e.g. the `read_run` result with `sample_accession=<biosample> AND library_strategy=WGS AND instrument_platform=ILLUMINA`), pick
   ~300 runs over >= 20 genera and >= 3 phyla, stream only the first N read pairs from the ENA FASTQ FTP (no full download).
   The label is the assembly's; the genome itself must be excluded from the reference (prerequisite 1), otherwise the test is a
   lookup. S-M.
2. **Coverage titration** 0.1x-50x by subsampling (first N reads, N from the genome size): accuracy and called fraction by rank
   versus depth; the minimum depth for genus and for species. S.
3. **Paired vs single and lane merging:** R1 only vs R1+R2 merged (the CLI merges; compare with concatenated input). S.
4. **Real vs simulated:** compare with `scripts/sim_reads.py` (substitution-only) at the same depths; add indels and
   quality-dependent errors. S.
5. **Cutoff generalisation:** the 0.99 reads cutoff was checked on one species pair; repeat leave-species-out for several
   genera with close species (Aspergillus section Fumigati, Fusarium oxysporum species complex, Candida/Nakaseomyces,
   Cryptococcus) to find the pairs it would name wrongly. M. This is the open risk of #32.

**Acceptance:** >= 5 genera from >= 3 labs, depth curves with the minimum depth for genus/species, paired-vs-single table.

## #43 RNA-seq

**Done.** 345 non-empty BFD RNA-seq sets at 100k/500k/2M reads (BFD's RNA files are partly 0-byte placeholders): median
containment .92-.95; genus precision .997 (called .93-.94), species called 63-68% at precision .966-.968 (`rna_reads_n*`).

**To add.**
1. **Stratify** by library type (poly-A, total RNA, strand-specific) and rRNA/host content: read the run metadata from the BFD
   download tables; if absent, infer from rRNA fraction (map the first 100k reads against the SILVA/rDNA hashes already in the
   index? simpler: count hashes shared with bacterial/rDNA decoys). S-M.
2. **The RNA-aware index (#33).** Most RNA k-mers come from transcribed regions (~50% of the genome in compact fungi, much less in
   intron-rich ones), so containment of the query in the *genome* is capped, which is why species calls plateau at ~65%. Build an
   index sketched from annotated transcripts/CDS (BFD has the annotations) and compare: genome index vs CDS index vs both. M-L; the
   largest unknown here.
3. **Contaminant reads:** add bacterial/host reads at 1-30% (from the decoy genomes) and measure the loss. S.
4. Spliced vs unspliced is not an issue for 21-mers inside exons; across junctions it costs ~20 k-mers per junction: measure on
   an intron-rich species (Cryptococcus) rather than assume. S.

**Acceptance:** accuracy vs depth by library type; genome vs CDS index comparison with a recommendation; contamination curve.

## #44 Long reads

**Done (2026-10-30): `reports/validation/long_reads.md`.** Simulated sweep (24 genomes) and 20 real runs; decision: a `high_error_reads` flag from the k-mer survival, no error-aware cutoff, no second index, `--min-abund auto` recommended at >= 10x.

**Problem.** A 21-mer survives an error rate e with probability (1-e)^21: at 5% 34%, at 10% 11%. But the containment *ANI*
(containment^(1/21)) tracks identity: a read set at 10% error looks like a 0.90 ANI relative, so it can reach a genus call at the
0.90 cutoff but never a species call at 0.9938/0.99, and the abundance filter (min count 2) throws away most true hashes at low coverage.

**Data.** Real runs with known species: ONT R10.4 / Kit14 (~1% error), older R9.4 (5-10%), PacBio HiFi (<0.5%) and CLR (10-15%); ~20
runs over 8+ species (NCBI SRA, `platform=OXFORD_NANOPORE|PACBIO_SMRT`; the species with chromosome-level assemblies from the same
BioSample give both the truth and the reference to exclude). **Simulation:** PBSIM3/Badread-style error sweep 0.5-15% with realistic
indel:substitution ratios on 10 genomes, depths 1x-30x.

**Questions and candidate changes (decide from the sweep):** error-aware species cutoff (estimate e from the abundance histogram
and shift the cutoff by the expected ANI loss); abundance filter that scales with coverage; a smaller k (15-17) or smaller
`scaled` for long reads only (needs a second index; k is fixed by the sketch). **Acceptance:** a table of accuracy vs error rate and
depth; a decision recorded for k / scaled / cutoff; `gmg run` flags `high_error_reads` when the histogram says so. **Effort:** L.
**Decision needed:** whether long reads are in scope for v0.1 (they are in the issue; I would put them after #42/#43).

## #45 Mixed samples and contaminated reads

**Depends on #31 (gather-style decomposition).** Without it the only possible statements are the ones already measured (the tool
reports the best single lineage; a foreign minority can win when the host's own reference is distant; no reliable mixture flag from
phylum-exclusive hashes). Design the experiment now so it runs the day #31 lands:

* Mixtures of 2 fungi at 50:50, 80:20, 90:10, 95:5, 99:1 (and 1:99 for the minor) for pairs at four distances (same species, same genus,
  same order, different phylum), assemblies and reads (5x-30x).
* Fungus + bacteria (E. coli, B. subtilis), fungus + human reads at 1-50%.
* Pure samples (negative: the detector must stay silent); cross-contamination at 0.1-1%.
* **Metrics:** dominant taxon correct (by rank); minor component detected (recall) as a function of its fraction; false mixture rate on
  pure inputs; abundance estimate error. Detection limit = the smallest fraction with recall >= 0.9 at false-rate <= 1%.

**#31 itself is an algorithm-design task** (greedy set cover of the query's hashes over the reference with abundance weighting, a
stopping rule from the error floor, per-component confidence). That is where I recommend a deeper model for the design review before
implementing. **Effort:** L (the experiment is S once the method exists).

## #46 Comparison with other tools

**Tools and how to run each fairly (same reference genomes, same truth sets):**

| tool | task | how | notes |
|---|---|---|---|
| skani (0.3.2 in the `bio` env) | nearest genome by ANI | `skani search` against a db of the BFD genomes (or its representatives) | the accuracy ceiling for assemblies; slower and larger index; label of the top hit |
| FastANI | same | bioconda | slow at 23k references; use the representatives set |
| mash screen / mash dist | containment of reads or genomes | `mash screen` against the BFD sketches (`results/mash_sketch` exists) | the closest in spirit; time/memory |
| sourmash gather (4.9.4) | decomposition, assemblies and reads | `sourmash gather` against a BFD sourmash database | also the reference for #31/#45; GmG sketches are sourmash-compatible |
| Kraken2 | k-mer LCA classifier | custom database from BFD with the NCBI taxids (now stored in the index metadata, `docs/taxonomy.md`) | needs a build (hours, tens of GB RAM); reads and contigs |
| ITS/rDNA classifier | marker-gene identification | extract ITS from the assemblies (ITSx) and classify with UNITE (vsearch/SINTAX) | taxonomy differs (UNITE vs NCBI): map by name; unavailable for reads without ITS |

**Sets:** post-freeze assemblies (#40), type-material (#39), Cocci reads (#42), RNA-seq (#43), stress inputs (#41), negative controls.
**Metrics:** accuracy by rank (same truth, same scoring, withheld = not called), wall time and peak memory per input, index size, and
failure modes: which inputs each tool gets wrong that the others get right. **Deliverable:** one comparison table and a short
"when to use which" note; a reproducible script per tool in `scripts/compare/`. **Effort:** L (the Kraken2 build and the ITS mapping
are the expensive parts; skani, mash and sourmash are S each). **Decision needed:** which of the six are worth the effort; I would do
skani, mash screen and sourmash gather first (they exist in the `bio` environment) and Kraken2 second.

## Order of work

1. Prerequisites 1-2 (S-M), then #39 (the truth set that sharpens every other number).
2. #40 D2, D3, D6 (S each: public data, no access issues) to reach six datasets; #41 remaining stress cases (S each).
3. #42 harvest and titration; #43 stratification; the cutoff-generalisation check for #32.
4. #46 with skani / mash / sourmash.
5. #43 RNA-aware index (#33), #44 long reads, #45 after #31.

## Decisions needed from you

* Where the lab assemblies (Coccidioides pangenome, MAT genomes, desert isolates) and the earlier mash IDs live (#40 D5).
* Whether a JGI Mycocosm export is available (#40 D4).
* Whether long reads (#44) and Kraken2/ITS comparisons (#46) are in scope for v0.1.
* Whether to build a small-input mode for mitochondrial/plasmid-sized inputs (#41), after measuring.
