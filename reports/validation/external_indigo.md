# External set: IndigoAg_Public endophyte/environmental assemblies (#40)

199 SPAdes assemblies of endophytic and environmental fungi (mostly Ascomycota: Sordariomycetes, Dothideomycetes, Leotiomycetes), none selected for BFD. Run: `scripts/hpcc/19_external_indigo.sbatch` (production index 2026.10.1 packed, balanced default), scored by `scripts/score_by_name.py` against the **provider's 2017 strain names** mapped to NCBI lineages through the taxdump (genus-first names only; species not scored). Per-genome tables stay on the HPCC.

Read this with care: the truth is a submitter's name, not curated taxonomy. Several "wrong" genus calls are legitimate modern splits (*Pestalotiopsis* -> *Pseudopestalotiopsis*, *Acremonium* -> *Sarocladium* / *Furcasterigmium*, *Lecanicillium* -> *Simplicillium* / *Akanthomyces*) and the provider genus is polyphyletic in places, so genus precision here is a lower bound; phylum and class are the more reliable rows.


status: {'ok': 138, 'low_confidence': 60, 'no_match': 1}; 25 without a usable genus name in the taxdump.

| rank | with truth | called | correct | precision | called share |
|---|---|---|---|---|---|
| phylum | 174 | 173 | 167 | 0.965 | 0.994 |
| class | 173 | 171 | 148 | 0.865 | 0.988 |
| order | 173 | 169 | 139 | 0.822 | 0.977 |
| family | 171 | 153 | 123 | 0.804 | 0.895 |
| genus | 174 | 116 | 79 | 0.681 | 0.667 |

Wrong phylum calls (truth, call): Ascomycota->Neocallimastigomycota x4, Ascomycota->Basidiomycota x2

Wrong genus calls: 37. Most frequent (provider genus -> call): Acremonium->Tetracladium, Lulworthia->Ramulariopsis, Lulworthia->Ramulariopsis, Lecanicillium->Simplicillium, Epichloe->Beauveria, Tiarosporella->Lizonia, Pestalotiopsis->Pseudopestalotiopsis, Epichloe->Xenoacremonium, Acremonium->Furcasterigmium, Acremonium->Sarocladium, Pestalotiopsis->Pseudopestalotiopsis, Magnaporthe->Colletotrichum, Lecanicillium->Akanthomyces, Pestalotiopsis->Pseudopestalotiopsis, Lophiostoma->Penicillium

## Reading

* **All six wrong phylum calls are `low_confidence`** with a top-hit ANI of 0.78-0.84 (four name Neocallimastigomycota, two Basidiomycota) for genera (*Lophodermium*, *Amphiporthe*, *Camarops*, *Hypocrea*, *Sarocladium*, *Coniochaeta*) that are thin or absent in BFD. This is the known behaviour of the phylum cutoff on divergent queries (#58); the status already says "do not trust". *Hypocrea* and *Sarocladium* belong to Hypocreales, so a 0.84 best hit in another phylum for them is odd and worth checking for contamination in those assemblies.
* 138 of 199 reach `ok` (a genus or species call); 60 are `low_confidence`; 1 `no_match`.
* The 25 inputs without a usable genus name are not scored.
