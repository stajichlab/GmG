# Selected strain sets from PopSeq (2026-11-13)

Four selections of public Illumina runs from PopSeq, made for popsketch population-structure pages. They are lists of runs, not
data: nothing has been downloaded or sketched yet.

| File | What | Runs | Countries | Studies |
|---|---|---|---|---|
| `saccharomyces_500.csv` | *S. cerevisiae* over how and where the strain was isolated | 500 | 67 | 141 |
| `afumigatus_pool150.csv` | *A. fumigatus*, a geographic pool to cut to 100 once clades are known | 150 | 33 | 75 |
| `coccidioides_immitis_50.csv` | *C. immitis* | 50 | 5 | 7 |
| `coccidioides_posadasii_50.csv` | *C. posadasii* | 50 | 7 | 12 |

Each `*.summary.tsv` gives, for every category, region, country, place and study, how many candidates there were and how many were
chosen. The columns of the CSVs are those of `popseq runs` plus `category`, `cultivation`, `keyword` (the word that decided the
category), `region`, `geo_loc_name` and `isolation_source` from the BioSample.

## How they were made

`make_popsets.sh` (run on the HPCC) does, in order: `popseq runs --species ... --layout PAIRED --min-depth 15` (PopSeq 0.2.0; ENA
Illumina paired FASTQ, not mutant, evolved, cross or non-WGS); `popset_attrs.py` joins the BioSample attributes PopSeq has harvested
into its cache (`isolation_source`, `env_*`, `geo_loc_name`, `host`, `strain` ...); `scripts/popset_select.py` (tests:
`scripts/test_popset_select.py`) chooses. The rules are in the docstring of `popset_select.py`; in short:

* one run per BioSample; estimated depth at least 25x (the coverage tier takes 25x), and among equal candidates the shallower run
  (it downloads less); at most 40 runs from one study (30 for *A. fumigatus*, 12 for each *Coccidioides* species);
* samples whose text says laboratory strain, evolution experiment, chemostat, FACS, spore, tetrad, plasmid, a cross, S288C and so on
  are left out: 18,067 of 69,276 *S. cerevisiae* candidates;
* *S. cerevisiae*: a category from keywords in the isolation text (quotas in the script: wine 70, beer 60, sake/rice 40, soy and
  other fermented food 40, bread and commercial 50, spirits and bioethanol 50, dairy 25, clinical 60, wild and nature 70, other
  known 20, no text 15; dairy had only 14, and what a small category cannot fill is shared by the others), then inside a category the
  run from the least represented country, continent and study; ties by a hash of the accession (seed 1), so the same input gives the
  same set;
* *A. fumigatus*: equal shares per continent, spread over countries and studies;
* *Coccidioides*: one species per file, spread over `geo_loc_name` (state or province) and studies.

## What these sets are not (read before using them)

* **The categories are a first pass by keywords.** Only 25% of *S. cerevisiae* BioSamples (17,884 of 71,311 candidates) have an
  isolation source (41% for *A. fumigatus*, 42% and 23% for *C. immitis* and *C. posadasii*). 43,545 candidates have no text at all
  and are `unannotated`; the 17 chosen are there to say so, not to represent them. `category` can be wrong: look at `keyword` and
  `isolation_source`. "Soy sauce" yeast is in `soy_other_ferment` only where a BioSample says soy (a few strain names do).
* **Species labels are NCBI's.** The earlier *Coccidioides* work found two label swaps and three non-*Coccidioides* samples in 559
  FASTQs, and the *C. posadasii* list here has samples from Michigan, Wisconsin, Minnesota and Germany (not endemic areas: clinical or
  mislabelled, not checked). Check the species gap in the sketches before using the labels.
* ***A. fumigatus* clades are not in PopSeq**, so "3 major clades" cannot be selected from metadata. This is a pool of 150 spread over
  continents to be sketched first; the final 100 should be chosen from the popsketch clusters (at least a third of the set in each
  clade), keeping the geographic spread. The pool is not clade-balanced.
* *Coccidioides* is almost all from the USA (immitis 45 of 50, posadasii 33 of 50): that is what the public data holds, the sets spread
  over states and studies, not over the world. There are only 86 *C. immitis* candidates at 25x, so 50 is more than half of them.
* Regions come from the run's country code; 5 *S. cerevisiae* and 23 *A. fumigatus* runs have none (`unknown_region`).
* The download is the full `.sra` or FASTQ of each run (the sketch is made at 25x, but the transfer is the whole run): the four
  sets hold 973, 440, 173 and 202 Gb of bases, about 0.4 bytes per base as an `.sra`, so about 0.4 TB, 0.18 TB and 0.07 TB each. No
  size ceiling is applied yet (PopSeq issue #48).

## Regenerating

    POPSEQ_ROOT=~/projects/PopSeq_prof_test PYTHON=~/projects/PopSeq/.venv/bin/python bash analysis/popsets/make_popsets.sh WORKDIR

reads the candidates and the biosample cache live, so a later PopSeq refresh can change the result; the CSVs here are the output of
2026-11-13. `scripts/popsra_fetch_chunks.slurm` takes a CSV of this kind as a manifest (it needs a `run_accession` column).

## Proposed to PopSeq (rule: metadata live there)

`isolation_source`, `env_broad_scale`, `env_local_scale`, `env_medium`, `lat_lon` and `geo_loc_name` are in PopSeq's biosample cache but not in
`popseq runs` or the panels; this folder reads the cache directly in `popset_attrs.py` as a stopgap. A curated `curated_cultivation`
column (spec 003 allows `curated_*` columns) would make the categories shared and correctable.
