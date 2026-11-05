#!/bin/bash
# rerun the Bd k-mer PCA without samples whose estimated genome is < 17 Mb (half the k-mer content), then recompare to the VCF PCA
set -euo pipefail
cd /bigdata/stajichlab/jstajich/projects/GmG_popgen/popgen_bd
EX=$(python3 - <<PY
import pandas as pd
q = pd.read_csv("out_s200/qc.tsv", sep="\t"); m = pd.read_csv("manifest.tsv", sep="\t").set_index("id")
b = q[(q.status == "ok") & (q.est_genome_mb < 17)]
import sys; print(f"{len(b)} samples < 17 Mb: lineage mix {m.organism.reindex(b.id).value_counts().to_dict()}; median cov {b.median_kmer_cov.median():.0f}, median min_abund {b.min_abund.median():.0f}; note mix {m.notes.reindex(b.id).value_counts().to_dict()}", file=sys.stderr)
print(",".join(b.id))
PY
)
python3 02_popdist.py sk200 manifest.tsv out_s200_qc --exclude "$EX" --cut 0.01 2>&1 | grep -v -i warn
KDIR=out_s200_qc python3 13b_bd_vcf_pca.py 2>&1 | grep -v -i warn
