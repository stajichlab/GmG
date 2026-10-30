#!/bin/bash -l
#SBATCH -J gmg_truth_prep
#SBATCH -p short
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/21_truth_prep.%j.log
# #39: truth lineage (pinned taxdump) for the RefSeq fungal assemblies, link to the index genomes, accession batches for download.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
REPO=${REPO:-$HOME/projects/GmG}
T=$GR/truthset; cd $T
REL=$GR/release/2026.10.2
TD=$GR/taxdump_20261029
python3 - <<'PY'
import json
with open("refseq_meta.tsv", "w") as o:
    o.write("id\ttaxid\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
    for line in open("refseq_fungi.jsonl"):
        j = json.loads(line)
        o.write(f"{j['accession']}\t{j['organism']['tax_id']}\t\t\t\t\t\t{j['organism'].get('organism_name','')}\n")
PY
python3 $REPO/scripts/normalize_taxonomy.py refseq_meta.tsv $TD refseq_lineage.tsv --report refseq_lineage.md
python3 $REPO/scripts/make_truth_set.py refseq_fungi.jsonl $REL/meta_full.tsv refseq_lineage.tsv truth
python3 $REPO/scripts/make_truth_sheets.py truth.all.tsv $REL/meta_full.tsv sheets
rm -rf batches; mkdir batches
tail -n +2 truth.all.tsv | cut -f1 | split -l 25 -d -a 3 - batches/batch_
ls batches | wc -l
