#!/bin/bash -l
#SBATCH -J gmg_thalia
#SBATCH -p short
#SBATCH -c 4
#SBATCH --mem=16G
#SBATCH -t 00:30:00
#SBATCH -o /bigdata/stajichlab/shared/projects/GmG_runs/logs/28_thalia.%j.log
# Do Conidiobolus-type fungi hit the Thalia democratica decoy (GCF_965202585.1) at ANI >= 0.90 (not_fungal)? The 10 ZyGoLife inputs
# whose top hit was that decoy, classified with and without the decoy (--exclude-ids), newest index.
set -uo pipefail
source $HOME/projects/GmG/scripts/hpcc/env.sh
REL=$GR/release/2026.10.3; [ -s $REL/bfd-2026.10.3-s10000.gmgi ] || REL=$GR/release/2026.10.2
INDEX=$(ls $REL/bfd-20*-s10000.gmgi)
R=$GR/reads/thalia; mkdir -p $R/in; cd $R
Z=/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Assembly/input
A=$GR/external/zygolife/in
for s in ARSEF_2338 ARSEF_6913 ARSEF_79 ARSEF_8428 NRRL_1229 NRRL_1231 NRRL_3011 NRRL_3624 NRRL_3643; do ln -sf $(readlink -f $A/$s.fa.gz) in/$s.fa.gz; done
echo GCF_965202585.1 > exclude.txt
echo "index $INDEX"
$GMG --threads 4 run in --index $INDEX -o with_decoy --no-report -q
$GMG --threads 4 run in --index $INDEX -o without_decoy --no-report -q --exclude-ids exclude.txt
# reads: the Conidiobolus antarcticus sample of the titration
head -1 $GR/reads/titration/sheet.tsv > sheet.tsv; grep "^Conidiobolus_antarcticus_ARSEF_691" $GR/reads/titration/sheet.tsv >> sheet.tsv
$GMG --threads 4 run --sample-sheet sheet.tsv --index $INDEX -o reads_with --no-report -q --max-reads 1500000
$GMG --threads 4 run --sample-sheet sheet.tsv --index $INDEX -o reads_without --no-report -q --max-reads 1500000 --exclude-ids exclude.txt
python3 - <<'PY'
import csv
for d in ("with_decoy", "without_decoy", "reads_with", "reads_without"):
    print("==", d)
    for r in csv.DictReader(open(f"{d}/results.tsv"), delimiter="\t"):
        print(r["query_id"][:34], r["status"], r["call_phylum"][:20], r["call_genus"][:16], r["ani"], r["top_hit_species"][:26], r["qc_flags"])
PY
