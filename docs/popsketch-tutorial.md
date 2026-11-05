# popsketch tutorial: audit a set of read files

This walks through one run on a simulated dataset that ships with the tool, so you can see what each output says before you point it at your own reads. The reference page is [popsketch.md](popsketch.md).

## 1. Make the example

```
popsketch demo -o demo/
```

This writes 32 small FASTQ files to `demo/reads/` and a label table to `demo/labels.tsv`. The reads are simulated from three groups of related genomes (9 samples each, groups about 1.5% divergent, members about 0.3% apart), plus five samples with a planted problem:

| sample | what is wrong |
|---|---|
| `lowcov0`, `lowcov1` | 0.7x k-mer coverage |
| `foreign` | reads of an unrelated genome |
| `half` | reads from only half of the genome |
| `mixed` | a second genome mixed in at 0.8x the host coverage |

and one wrong label: `g0s0` is labelled `group1` in `labels.tsv` but its reads are group 0. These are the problems a real set can hide; the example lets you check that you can see them.

## 2. Run it

```
popsketch run demo/reads -o demo/out --labels demo/labels.tsv --scaled 5 --max-bases 0
```

`--scaled 5` and `--max-bases 0` are for these tiny genomes (120 kb); for a real fungal genome leave both at their defaults (scaled 200, 1 Gb of bases per sample). The run takes a few seconds and ends with a line such as `32 samples, 5 flagged, 27 in the analysis`. Open `demo/out/report.html`.

## 3. Read the audit first

The report opens with the audit, and `qc.tsv` has the same table for every sample. The flagged rows of the example, with the reasons the tool gives:

```
foreign   excluded  unrelated_to_cohort  nearest neighbour ANI 0.0000 vs cohort median 0.9970 (z -199.4): shares almost nothing with the other samples
half      excluded  size_low             estimated size 0.06 Mb vs cohort median 0.12 Mb (z -16.5); smaller than the cohort: content missing
lowcov0   excluded  low_coverage         k-mer coverage 0.7 is below 4.0
lowcov1   excluded  low_coverage         k-mer coverage 0.6 is below 4.0
mixed     excluded  size_high            estimated size 0.24 Mb vs cohort median 0.12 Mb (z +33.0); larger than the cohort: a second genome or extra content
```

Exactly the five planted samples are flagged and nothing else. Every flag is relative to the rest of the set, so nothing here needs a reference. By default flagged samples are left out of the distances, PCA, clusters and tree (`--flagged exclude`); `--flagged mark` keeps them in with their own marker, and `--flagged decontaminate` cleans samples flagged only for contamination and puts them back.

What the audit cannot see is stated in the report: light contamination (about 1 to 7% of the reads) is invisible without `--index`.

## 4. Check the labels

`label_check.tsv` lists samples whose ANI matches another label better than their own:

```
group   g0s0   group1   own_median_ani 0.9677   better_matching_label group0   its_median_ani 0.9969   nearest_sample g0s1 (group0, 0.9970)
```

`g0s0` is closer to group 0 (0.9969) than to the group it is labelled with (0.9677), and its nearest neighbour is in group 0. It is the swapped label. This is evidence to check the sample's metadata, not proof of an error. The second label column in the example, `source`, has blanks (unlabelled) and is not a genetic grouping; the tool does not check label groups that do not form coherent clusters in the sketches, so it reports no conflicts for it and says why.

## 5. Then the structure

For the 27 samples that passed: the PCA (`pca.tsv`, plotted in the report), clusters, and a tree (`tree.nwk`). `popsketch.json` records how the clusters were chosen:

```
verdict clear:3, the jump is 13.3 times the noise floor, null probability 0.0043
noise floor 0.00214 (measured from split halves of the samples)
```

The cut into 3 clusters was made because the largest jump in the clustering is 13.3 times the technical noise floor (measured by comparing two halves of the same reads). `clusters.tsv` gives the cluster of each sample; against the group labels the adjusted Rand index is 0.887 (the swapped sample is the one disagreement). If no jump had cleared 3 noise floors, the report would say "no single clear cut" and list candidates for `--k` instead of forcing one.

The last section of the report states what these data can and cannot say. In this example no pair is closer than the noise floor, so within-group structure is weakly supported; in a near-clonal set most pairs fall inside the noise and the report says to use a SNP-based method for strain-level relationships.

## 6. On your own data

- Point `popsketch run` at a folder of FASTQ files (R1/R2 and lanes of one sample are merged by their names), or give a sample sheet (`--sample-sheet`: tab-separated, columns `id` and `path`, several rows with one id are one sample).
- Make a label table: tab-separated, first column the sample id, then one column per label set (species, lineage, source...); leave a cell blank when unknown. The first label column colours the plots.
- Add `--index bfd.gmgi` (a GmG reference index) to name the organism of each sample and to catch mixtures at 1 to 7% that the abundance spectrum cannot see.
- If you already have a SNP tree or PCA of the same samples, `--compare-tree` and `--compare-pcs` report how well the sketches agree with them, and whether each cluster is a clade of your tree.
- `--report-format md` or `all` writes `report.md` with each figure as an SVG file.
- Sets of up to 1,500 samples are supported; above `--max-samples` the run stops with a message instead of running for hours.
