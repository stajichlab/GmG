# Stress tests: degraded known genomes

194 inputs from 24 chromosome-level assemblies (genus present in BFD, not in the index). Truth = the host lineage. Cells: called/total (correct among called). Species is not scored: 23 of the 24 host assemblies have no species name in NCBI, so the species columns only show that a call was made.

| variant | param | n | phylum called (correct) | genus called (correct) | species called (correct) | median hashes | flags seen |
|---|---|---|---|---|---|---|---|
| complete | 10 | 24 | 23/24 (23) | 19/24 (19) | 16/24 (1) | 453 | short_genome×14 |
| complete | 25 | 24 | 24/24 (24) | 19/24 (19) | 16/24 (1) | 1024 |  |
| complete | 50 | 24 | 24/24 (24) | 19/24 (19) | 16/24 (1) | 1947 |  |
| complete | 75 | 24 | 24/24 (24) | 19/24 (19) | 16/24 (1) | 2894 |  |
| frag5k | 5000 | 24 | 24/24 (24) | 19/24 (19) | 16/24 (1) | 3790 | fragmented×24 |
| contam | 5 | 24 | 24/24 (21) | 19/24 (19) | 15/24 (1) | 3994 |  |
| contam | 20 | 24 | 24/24 (19) | 24/24 (19) | 0/24 (0) | 4698 |  |
| contam | 50 | 24 | 24/24 (17) | 23/24 (16) | 0/24 (0) | 7041 |  |

## Negative controls

* `neg_random_acgt`: status `no_match`, phylum call `-`, flags ``
* `neg_repeat`: status `no_match`, phylum call `-`, flags `few_hashes`

A tentative genus is not a call; `call_*` columns only. Wrong calls are those where the called taxon differs from the host lineage; for the contamination series the host is the majority genome.


## Reading

* **Completeness is not the limit.** Down to 10% of a genome (median 453 hashes) the phylum is right for 23/24 and genus calls are the same 19/19 correct as at 75%; the 5 uncalled hosts are divergent genera that are already `low_confidence` at full length. The `short_genome` flag fires at 10% for 14 of 24.
* **Fragmentation does not matter for the sketch.** 5 kb contigs give the same calls as the intact genome; the `fragmented` flag fires for all 24 (the QC says so; the call is unaffected).
* **Contamination from another phylum is the weak point.** Phylum correct 21/24 with 5% foreign DNA, 19/24 at 20%, 17/24 at 50% (genus wrong for 5 and 8 of those called at 20% and 50%). The three failures at 5% are Basidiomycota hosts whose nearest BFD relative is already distant (uncontaminated ANI 0.77-0.85, no genus call): a 5% Ascomycota donor that has a close relative in BFD outcompetes them and the top hit flips to Ascomycota (ANI 0.87-0.88). A single-lineage call cannot say "this is a mixture"; the host in these cases is the minority-hash source.
* **Negative controls pass**: uniform random DNA and a tandem repeat give `no_match` with no phylum call (the repeat also flags `few_hashes`).
* Not covered here: bacterial/human contamination (no such genomes were fetched), chimeras of two fungi from the same phylum, hybrids, mitochondrion/plasmid-only inputs.
