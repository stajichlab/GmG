# The non-trivial species_match verdicts, looked at (2026-11-18)

Read evidence: the first 3M read pairs of each run (`10_classify_changed.sbatch`: minimap2 -x sr to the filed species' reference, and `gmg guess`) plus PopSeq catalog metadata. A verdict says what the reads are; the reasons below are checked only where the evidence is named. Nothing has been changed in the catalog.

| run | filed | verdict | what the reads are | evidence | reading |
|---|---|---|---|---|---|
| SRR30534096 | *C. neoformans* | mismatch | *A. fumigatus* | 27,014 hashes (33%); 0.2% in the filed genus; 0.8% map to H99 | mislabelled |
| SRR34066268, SRR34066269 | *C. neoformans* | genus_match | *C. deneoformans* | 86 to 93% of hashes in *Cryptococcus*, species call *deneoformans* | species mislabel; strains CBS7816 and CBS7748, same study (PRJNA1279655, "Cryptococcus neoformans genome sequencing") as the runs in `../RESULTS_hybrids_and_155.md`; reads not mapped to JEC21 yet |
| ERR12084585, SRR31166608, SRR31166716 | *C. beticola* | mixed | *Cercospora* plus *Penicillium* or *Akanthomyces* | 32 to 50% of hashes in *Cercospora*, 47 to 58% in the other genus; 5 to 34% of reads map to *C. beticola* | mixed samples (probably leaf material), not mislabels; not confirmed from sample metadata |
| ERR12084630 | *C. beticola* | host_dominated | not called | purity 14%; 4.8% of reads map to *C. beticola* | mixed or plant material, same kind as the `mixed` runs |
| SRR31166705 | *C. beticola* | no_fungal_call_suspect | not called | purity 50%, no call | same kind as the mixed runs; not `host_dominated` because the filed genus holds more than `--mixed-purity` without a call: the index should have called it; open |
| ERR1976666 to ERR1976671 (PRJEB20844) | *P. oryzae* | host_dominated | wheat tissue with some *P. oryzae* | study title "First report of blast disease in wheat caused by *Magnaporthe oryzae* ..."; ERR1976668/670/671 call *P. oryzae* (560 to 1,389 hashes) but only 10%, 16%, 5% of reads map; ERR1976667/669 only weak hits (about 230 to 290 hashes) to rust-like genomes | infected-plant samples: the organism is right, the read set is mostly host. ERR1976666 mapping not available (reference download failed; retried) |
| SRR26462617 | *C. neoformans* | no_fungal_call_suspect | not fungal in this index | 0.03% map to the filed species' reference; purity 0.1% | probably non-fungal (not identified) |
| ERR13824414 | *Z. tritici* | host_dominated | not fungal in this index | 0.5% map to the filed species' reference; purity 1.4% | probably non-fungal, but at 1.4% the rule calls it host-dominated; the status does not say what dominates (not identified) |
| SRR10233430 | *C. neoformans* | host_dominated | not examined | purity 5%, best hit 1,001 hashes; the study (PRJNA575911) is titled "*Cryptococcus gattii* VGI Raw sequence reads" and the strain field says "Fungi" | the study title says *C. gattii*, the catalog species says *C. neoformans*: possibly mislabelled in the other direction, not checked against the reads |

## What follows for the check
- **Purity separates the cases that matter.** A mislabel has almost none of the filed genus (SRR30534096: 0.2%); mixtures hold a third to a half; infected-plant samples hold a few percent with a weak call. Only the first should be excluded by default, as agreed.
- **Plant-tissue runs are a class of their own, and `host_dominated` names it.** The six *P. oryzae* runs are filed correctly; the problem is that they are not isolate genomes. They are `host_dominated` (owner's decision 2026-11-18: excluded by default).
- Not examined against the reads: SRR10233430; the *C. deneoformans* pair against JEC21; the sample metadata of the *C. beticola* study.
